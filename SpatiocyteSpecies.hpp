//
// E-Cell is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public
// License as published by the Free Software Foundation; either
// version 2 of the License, or (at your option) any later version.
// 
// E-Cell is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public
// License along with E-Cell -- see the file COPYING.
// If not, write to the Free Software Foundation, Inc.,
// 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
// 
//END_HEADER
//
// written by Satya Arjunan <satya.arjunan@gmail.com>
// E-Cell Project, Institute for Advanced Biosciences, Keio University.
//


#ifndef __SpatiocyteSpecies_hpp
#define __SpatiocyteSpecies_hpp

#include <sstream>
#include <Variable.hpp>
#include <gsl/gsl_randist.h>
#include "SpatiocyteCommon.hpp"
#include "SpatiocyteStepper.hpp"
#include "SpatiocyteProcessInterface.hpp"
#include "DiffusionInfluencedReactionProcessInterface.hpp"
#include "MoleculePopulateProcessInterface.hpp"

// The size of Coord must be 128 bytes to avoid cacheline splits
// The Core 2 has 64-byte cacheline
static double getDistance(Point* aSourcePoint, Point* aDestPoint)
{
  return sqrt(pow(aDestPoint->x-aSourcePoint->x, 2)+
              pow(aDestPoint->y-aSourcePoint->y, 2)+
              pow(aDestPoint->z-aSourcePoint->z, 2));
}

String int2str(int anInt)
{
  std::stringstream aStream;
  aStream << anInt;
  return aStream.str();
}

/*
 * isVacant: vacant species definition: a species on which other species can
 * diffuse on or occupy. One major optimization in speed for vacant species is
 * that theCoords list is not updated when other species diffuse on it. 
 * There are four possible types of vacant species:
 * 1. !isDiffusiveVacant && !isReactiveVacant: a vacant species
 *    whose theCoords list and theCoordSize are never updated. This is
 *    the most basic version. This vacant species never diffuses, or reacts
 *    using SNRP. 
 * 2. !isDiffusiveVacant && isReactiveVacant: a vacant species which is a
 *    substrate of a SNRP reaction. In this case theCoordSize is updated
 *    before the step interval of SNRP is calculated, and theCoords list is
 *    updated before the SNRP reaction (fire) is executed to get a valid
 *    molecule list. Set by SNRP.
 * 3. isDiffusiveVacant && !isReactiveVacant: a vacant species which also
 *    diffuses. In this case, theCoords list and theCoordsSize are
 *    updated just before it is diffused. Set by DiffusionProcess.   
 * 4. isDiffusiveVacant && isReactiveVacant: a vacant species which reacts
 *    using SNRP and also diffuses.
 * 5. isCompVacant: the VACANT species declared for each compartment. It can
 *    also be isDiffusiveVacant and isReactiveVacant. Set during compartment
 *    registration. It also persistently stores all the compartment voxels.
 *    Referred to as theVacantSpecies. For isCompVacant, initially there
 *    are no molecules in its list. All voxels are stored in theCompCoords. Only
 *    if it is updated before being called by VisualizationLogProcess, SNRP or
 *    DiffusionProcess, it will be theCoords will be populated with the
 *    CompCoords.
 * 6. isVacant {isCompVacant; isDiffusiveVacant; isReactiveVacant): the
 *    general name used to identify either isCompVacant, isDiffusiveVacant or
 *    isReactiveVacant to reduce comparison operations.
 */

class Species
{
public:
  Species(SpatiocyteStepper* aStepper, Variable* aVariable, int anID, 
          int anInitCoordSize, const gsl_rng* aRng, double voxelRadius,
          std::vector<Voxel>& aLattice):
    isCentered(false),
    isCompVacant(false),
    isDiffusing(false),
    isDiffusiveVacant(false),
    isFixedAdjoins(false),
    isGaussianPopulation(false),
    isInContact(false),
    isOffLattice(false),
    isPolymer(false),
    isReactiveVacant(false),
    isSubunitInitialized(false),
    isVacant(false),
    theID(anID),
    theCollision(0),
    theInitCoordSize(anInitCoordSize),
    theCoordSize(0),
    D(0),
    theDiffusionInterval(libecs::INF),
    theWalkProbability(1),
    theRadius(voxelRadius),
    theRng(aRng),
    thePopulateProcess(NULL),
    theStepper(aStepper),
    theVariable(aVariable),
    theCompCoords(&theCoords),
    theLattice(aLattice) {}
  ~Species() {}
  void initialize(int speciesSize, int anAdjoiningCoordSize)
    {
      theAdjoiningCoordSize = anAdjoiningCoordSize;
      theReactionProbabilities.resize(speciesSize);
      theDiffusionInfluencedReactions.resize(speciesSize);
      theFinalizeReactions.resize(speciesSize);
      for(int i(0); i != speciesSize; ++ i)
        {
          theDiffusionInfluencedReactions[i] = NULL;
          theReactionProbabilities[i] = 0;
          theFinalizeReactions[i] = false;
        }
      if(theComp)
        {
          setVacantSpecies(theComp->vacantSpecies);
        }
    }
  void setDiffusionInfluencedReaction(
                                    DiffusionInfluencedReactionProcessInterface*
                                      aReaction, int anID, double aProbability)
    {
      theDiffusionInfluencedReactions[anID] = aReaction;
      theReactionProbabilities[anID] = aProbability;
    }
  void setDiffusionInfluencedReactantPair(Species* aSpecies)
    {
      theDiffusionInfluencedReactantPairs.push_back(aSpecies);
    }
  void setPopulateProcess(MoleculePopulateProcessInterface* aProcess,
                          double aDist)
    {
      if(aDist)
        {
          isGaussianPopulation = true;
        }
      thePopulateProcess = aProcess;
    }
  bool getIsGaussianPopulation()
    {
      return isGaussianPopulation;
    }
  int getPopulatePriority()
    {
      return thePopulateProcess->getPriority();
    }
  void populateCompGaussian()
    {
      if(thePopulateProcess)
        {
          thePopulateProcess->populateGaussian(this);
        }
      else if(theCoordSize)
        {
          std::cout << "Warning: Species " <<
            theVariable->getFullID().asString() <<
            " not populated." << std::endl;
        }
    }
  bool getIsPopulateSpecies()
    {
      return (thePopulateProcess != NULL);
    }
  void populateCompUniform(unsigned int voxelIDs[], unsigned int* aCount)
    {
      if(thePopulateProcess)
        {
          thePopulateProcess->populateUniformDense(this, voxelIDs, aCount);
        }
      else if(theCoordSize)
        {
          std::cout << "Species:" << theVariable->getFullID().asString() <<
            " not CoordPopulated." << std::endl;
        }
    }
  void populateCompUniformSparse()
    {
      if(thePopulateProcess)
        {
          thePopulateProcess->populateUniformSparse(this);
        }
      else if(theCoordSize)
        {
          std::cout << "Species:" << theVariable->getFullID().asString() <<
            " not CoordPopulated." << std::endl;
        }
    }
  void populateUniformOnDiffusiveVacant()
    {
      if(thePopulateProcess)
        {
          thePopulateProcess->populateUniformOnDiffusiveVacant(this);
        }
      else if(theCoordSize)
        {
          std::cout << "Species:" << theVariable->getFullID().asString() <<
            " not CoordPopulated." << std::endl;
        }
    }
  Variable* getVariable() const
    {
      return theVariable;
    }
  std::vector<unsigned int> getSourceCoords()
    {
      std::vector<unsigned int> aCoords;
      for(unsigned int i(0); i != theCoordSize; ++i)
        {
          std::vector<unsigned int>& 
            aSourceCoords(theLattice[theCoords[i]].subunit->sourceCoords);
          for(unsigned int j(0); j != aSourceCoords.size(); ++j)
            {
              if(aSourceCoords[j] != theNullCoord)
                {
                  aCoords.push_back(aSourceCoords[j]);
                }
            }
        }
      return aCoords;
    }
  std::vector<unsigned int> getTargetCoords()
    {
      std::vector<unsigned int> aCoords;
      for(unsigned int i(0); i != theCoordSize; ++i)
        {
          std::vector<unsigned int>& 
            aTargetCoords(theLattice[theCoords[i]].subunit->targetCoords);
          for(unsigned int j(0); j != aTargetCoords.size(); ++j)
            {
              if(aTargetCoords[j] != theNullCoord)
                {
                  aCoords.push_back(aTargetCoords[j]);
                }
            }
        }
      return aCoords;
    }
  std::vector<unsigned int> getSharedCoords()
    {
      std::vector<unsigned int> aCoords;
      for(unsigned int i(0); i != theCoordSize; ++i)
        {
          std::vector<unsigned int>& 
            aSharedLipids(theLattice[theCoords[i]].subunit->sharedLipids);
          for(unsigned int j(0); j != aSharedLipids.size(); ++j)
            {
              if(aSharedLipids[j] != theNullCoord)
                {
                  aCoords.push_back(aSharedLipids[j]);
                }
            }
        }
      return aCoords;
    }
  std::vector<unsigned int>& getCoords()
    {
      return theCoords;
    }
  unsigned int size() const
    {
      return theCoordSize;
    }
  unsigned int getCoord(int anIndex)
    {
      return theCoords[anIndex];
    }
  Point getPoint(int anIndex)
    {
      if(isOffLattice)
        {
          return *theLattice[theCoords[anIndex]].point;
        }
      else if(isPolymer)
        {
          return theLattice[theCoords[anIndex]].subunit->subunitPoint;
        }
      return theStepper->coord2point(theCoords[anIndex]);
    }
  unsigned short getID() const
    {
      return theID;
    }
  double getMeanSquaredDisplacement()
    {
      if(!theCoordSize)
        {
          return 0;
        }
      double aDisplacement(0);
      for(unsigned int i(0); i < theCoordSize; ++i)
        {
          Point aCurrentPoint(theStepper->getPeriodicPoint(theCoords[i],
                                                 theDimension,
                                                 &theCoordOrigins[i]));
          double aDistance(getDistance(&theCoordOrigins[i].point,
                                       &aCurrentPoint));
          aDisplacement += aDistance*aDistance;
        }
      return
        aDisplacement*pow(theRadius*2, 2)/theCoordSize;
    }
  void setCollision(unsigned int aCollision)
    {
      theCollision = aCollision;
    }
  void setIsSubunitInitialized()
    {
      isSubunitInitialized = true;
    }
  void setIsCompVacant()
    {
      isCompVacant = true;
      isVacant = true;
      setVacantSpecies(this);
    }
  void setIsInContact()
    {
      isInContact = true;
    }
  void setIsCentered()
    {
      isCentered = true;
    }
  void setIsPopulated()
    {
      theInitCoordSize = theCoordSize;
      getVariable()->setValue(theCoordSize);
    }
  void finalizeSpecies()
    {
      if(theCollision)
        {
          collisionCnts.resize(theCoordSize);
          for(std::vector<unsigned int>::iterator 
              i(collisionCnts.begin()); i != collisionCnts.end(); ++i)
            {
              *i = 0;
            }
        }
      //need to shuffle molecules of the compVacant species if it has
      //diffusing vacant species to avoid bias when random walking:
      if(isCompVacant)
        {
          for(unsigned int i(0); i != theComp->species.size(); ++i)
            {
              if(theComp->species[i]->getIsDiffusiveVacant())
                {
                  shuffleCoords();
                  break;
                }
            }
        }
    }
  void shuffleCoords()
    {
      std::random_shuffle(theCoords.begin(), theCoords.end());
    }
  unsigned int getCollisionCnt(unsigned int anIndex)
    {
      return collisionCnts[anIndex];
    }
  unsigned int getCollision() const
    {
      return theCollision;
    }
  void setIsDiffusiveVacant()
    {
      isDiffusiveVacant = true;
      isVacant = true;
    }
  void setIsReactiveVacant()
    {
      isReactiveVacant = true;
      isVacant = true;
    }
  void setIsOffLattice()
    {
      isOffLattice = true;
      isFixedAdjoins = false;
      for(unsigned int i(0); i != theComp->species.size(); ++i)
        {
          theComp->species[i]->setIsFixedAdjoins(false);
        }
    }
  void setIsFixedAdjoins(bool state)
    {
      isFixedAdjoins = state;
    }
  void setIsPolymer(std::vector<double> bendAngles, int aDirectionality)
    {
      theBendAngles.resize(0);
      thePolymerDirectionality = aDirectionality;
      isPolymer = true;
      for(std::vector<double>::const_iterator i(bendAngles.begin()); 
          i != bendAngles.end(); ++i)
        {
          theBendAngles.push_back(*i);
          if(thePolymerDirectionality != 0)
            {
              theBendAngles.push_back(*i+M_PI);
            }
        }
    }
  void setDiffusionCoefficient(double aCoefficient)
    {
      D = aCoefficient;
      if(D > 0)
        {
          isDiffusing = true;
        }
    }
  double getDiffusionCoefficient() const
    {
      return D;
    }
  double getWalkProbability() const
    {
      return theWalkProbability;
    }
  bool getIsPolymer() const
    {
      return isPolymer;
    }
  bool getIsOffLattice()
    {
      return isOffLattice;
    }
  bool getIsSubunitInitialized() const
    {
      return isSubunitInitialized;
    }
  bool getIsDiffusing() const
    {
      return isDiffusing;
    }
  bool getIsCompVacant() const
    {
      return isCompVacant;
    }
  bool getIsVacant() const
    {
      return isVacant;
    }
  bool getIsDiffusiveVacant()
    {
      return isDiffusiveVacant;
    }
  bool getIsReactiveVacant()
    {
      return isReactiveVacant;
    }
  bool getIsLipid() const
    {
      return (isCompVacant && theComp->dimension == 2);
    }
  bool getIsInContact() const
    {
      return isInContact;
    }
  bool getIsCentered() const
    {
      return isCentered;
    }
  bool getIsPopulated() const
    {
      return theCoordSize == theInitCoordSize;
    }
  double getDiffusionInterval() const
    {
      return theDiffusionInterval;
    }
  unsigned int getDimension()
    {
      return theDimension;
    }
  void setDimension(unsigned int aDimension)
    {
      theDimension = aDimension;
      if(theDimension == 3)
        {
          isFixedAdjoins = true;
        }
    }
  void resetFinalizeReactions()
    {
      for(unsigned int i(0); i != theFinalizeReactions.size(); ++i)
        {
          theFinalizeReactions[i] = false;
        }
    }
  void finalizeReactions()
    {
      for(unsigned int i(0); i != theFinalizeReactions.size(); ++i)
        {
          if(theFinalizeReactions[i])
            {
              theDiffusionInfluencedReactions[i]->finalizeReaction();
            }
        }
    }
  void addCollision(unsigned int aCoord)
    {
      for(unsigned int i(0); i < theCoordSize; ++i)
        {
          if(aCoord == theCoords[i])
            {
              ++collisionCnts[i];
              return;
            }
        }
      std::cout << "error in species add collision" << std::endl;
    }
  void collide()
    {
      for(unsigned int i(0); i < theCoordSize; ++i)
        {
          Voxel& source(theLattice[theCoords[i]]);
          int size;
          if(isFixedAdjoins)
            {
              size = theAdjoiningCoordSize;
            }
          else
            {
              size = source.diffuseSize;
            }
          Voxel& target(theLattice[source.adjoiningCoords[
                        gsl_rng_uniform_int(theRng, size)]]);
          if(target.id == theVacantID)
            {
              if(theWalkProbability == 1 ||
                 gsl_rng_uniform(theRng) < theWalkProbability)
                {
                  target.id = theID;
                  source.id = theVacantID;
                  theCoords[i] = target.coord;
                }
            }
          else if(theDiffusionInfluencedReactions[target.id])
            {
              //If it meets the reaction probability:
              if(gsl_rng_uniform(theRng) < theReactionProbabilities[target.id])
                { 
                  Species* targetSpecies(theStepper->id2species(target.id));
                  ++collisionCnts[i];
                  targetSpecies->addCollision(target.coord);
                  if(theCollision == 2)
                    {
                      DiffusionInfluencedReactionProcessInterface* aReaction(
                                 theDiffusionInfluencedReactions[target.id]);
                      if(aReaction->react(source.coord, target.coord))
                        {
                          //Soft remove the source molecule, i.e.,
                          //keep the id intact:
                          theCoords[i--] = theCoords[--theCoordSize];
                          theVariable->setValue(theCoordSize);
                          //Soft remove the target molecule:
                          targetSpecies->softRemoveCoord(target.coord);
                          theFinalizeReactions[targetSpecies->getID()] = true;
                        }
                    }
                }
            }
        }
    }
  void walk()
    {
      for(unsigned int i(0); i < theCoordSize; ++i)
        {
          Voxel& source(theLattice[theCoords[i]]);
          int size;
          if(isFixedAdjoins)
            {
              size = theAdjoiningCoordSize;
            }
          else
            {
              size = source.diffuseSize;
            }
          Voxel& target(theLattice[source.adjoiningCoords[
                        gsl_rng_uniform_int(theRng, size)]]);
          if(target.id == theVacantID)
            {
              if(theWalkProbability == 1 ||
                 gsl_rng_uniform(theRng) < theWalkProbability)
                {
                  target.id = theID;
                  source.id = theVacantID;
                  theCoords[i] = target.coord;
                }
            }
          else if(theDiffusionInfluencedReactions[target.id])
            {
              //If it meets the reaction probability:
              if(gsl_rng_uniform(theRng) < theReactionProbabilities[target.id])
                { 
                  Species* targetSpecies(theStepper->id2species(target.id));
                  DiffusionInfluencedReactionProcessInterface* aReaction(
                             theDiffusionInfluencedReactions[target.id]);
                  if(aReaction->react(source.coord, target.coord))
                    {
                      //Soft remove the source molecule, i.e.,
                      //keep the id intact:
                      theCoords[i--] = theCoords[--theCoordSize];
                      theVariable->setValue(theCoordSize);
                      //Soft remove the target molecule:
                      targetSpecies->softRemoveCoord(target.coord);
                      theFinalizeReactions[targetSpecies->getID()] = true;
                    }
                }
            }
        }
    }
  void walkVacant()
    {
      updateVacantCoords();
      for(unsigned int i(0); i < theCoordSize; ++i)
        {
          Voxel& source(theLattice[theCoords[i]]);
          int size;
          if(isFixedAdjoins)
            {
              size = theAdjoiningCoordSize;
            }
          else
            {
              size = source.diffuseSize;
            }
          Voxel& target(theLattice[source.adjoiningCoords[
                        gsl_rng_uniform_int(theRng, size)]]);
          if(target.id == theVacantID)
            {
              if(theWalkProbability == 1 ||
                 gsl_rng_uniform(theRng) < theWalkProbability)
                {
                  target.id = theID;
                  source.id = theVacantID;
                }
            }
          /*
          else if(theDiffusionInfluencedReactions[target.id] != NULL)
            {
              //If it meets the reaction probability:
              if(gsl_rng_uniform(theRng) <
                 theReactionProbabilities[target.id])
                { 
                  Species* targetSpecies(theStepper->id2species(target.id));
                  DiffusionInfluencedReactionProcessInterface* aReaction(
                             theDiffusionInfluencedReactions[target.id]);
                  if(aReaction->react(source.coord, target.coord))
                    {
                      //Soft remove the target molecule:
                      targetSpecies->softRemoveCoord(target.coord);
                      theFinalizeReactions[targetSpecies->getID()] = true;
                    }
                }
            }
            */
        }
    }
  void setComp(Comp* aComp)
    {
      theComp = aComp;
    }
  Comp* getComp()
    {
      return theComp;
    }
  void setVariable(Variable* aVariable)
    {
      theVariable = aVariable;
    }
  void removeSurfaces()
    {
      int newCoordSize(0);
      for(unsigned int i(0); i < theCoordSize; ++i) 
        {
          unsigned int aCoord(theCoords[i]);
          if(theStepper->isRemovableEdgeCoord(aCoord, theComp))
            {
              Comp* aSuperComp(
                 theStepper->system2Comp(theComp->system->getSuperSystem())); 
              aSuperComp->vacantSpecies->addCompCoord(aCoord);
            }
          else 
            { 
              theCoords[newCoordSize] = aCoord; 
              ++newCoordSize; 
            }
        }
      theCoordSize = newCoordSize;
      theVariable->setValue(theCoordSize);
    }
  void removePeriodicEdgeCoords()
    {
      int newCoordSize(0);
      for(unsigned int i(0); i < theCoordSize; ++i) 
        {
          unsigned int aCoord(theCoords[i]);
          if(theStepper->isPeriodicEdgeCoord(aCoord, theComp))
            {
              theLattice[aCoord].id = theLattice[theNullCoord].id;
            }
          else 
            { 
              theCoords[newCoordSize] = aCoord; 
              ++newCoordSize; 
            }
        }
      theCoordSize = newCoordSize;
      theVariable->setValue(theCoordSize);
    }
  void updateSpecies()
    {
      if(isCompVacant && (isDiffusiveVacant || isReactiveVacant))
        {
          theCompCoords = new std::vector<unsigned int>;
          for(unsigned int i(0); i != theCoordSize; ++i)
            { 
              theCompCoords->push_back(theCoords[i]);
            }
        }
    }
  //If it isReactiveVacant it will only be called by SNRP when it is substrate
  //If it isDiffusiveVacant it will only be called by DiffusionProcess before
  //being diffused. So we need to only check if it isVacant:
  void updateCoords()
    {
      if(isDiffusiveVacant || isReactiveVacant)
        {
          updateVacantCoords();
        }
    }
  //If it isReactiveVacant it will only be called by SNRP when it is substrate:
  void updateCoordSize()
    {
      if(isDiffusiveVacant || isReactiveVacant)
        {
          updateVacantCoordSize();
        }
    }
  //Even if it is a isCompVacant, this method will be called by
  //VisualizationLogProcess, or SNRP if it is Reactive or DiffusionProcess
  //if it is Diffusive:
  void updateVacantCoords()
    {
      theCoordSize = 0;
      int aSize(theVacantSpecies->compCoordSize());
      for(int i(0); i != aSize; ++i)
        { 
          unsigned int aCoord(theVacantSpecies->getCompCoord(i));
          if(theLattice[aCoord].id == theID)
            {
              ++theCoordSize;
              if(theCoordSize > theCoords.size())
                {
                  theCoords.push_back(aCoord);
                }
              else
                {
                  theCoords[theCoordSize-1] = aCoord;
                }
            }
        }
      theVariable->setValue(theCoordSize);
    }
  void updateVacantCoordSize()
    {
      theCoordSize = 0;
      int aSize(theVacantSpecies->compCoordSize());
      for(int i(0); i != aSize; ++i)
        { 
          unsigned int aCoord(theVacantSpecies->getCompCoord(i));
          if(theLattice[aCoord].id == theID)
            {
              ++theCoordSize;
            }
        }
      if(theCoordSize > theCoords.size())
        {
          theCoords.resize(theCoordSize);
        }
      theVariable->setValue(theCoordSize);
    }
  void addCoord(unsigned int aCoord)
    {
      theLattice[aCoord].id = theID;
      if(!isVacant)
        {
          ++theCoordSize;
          if(theCoordSize > theCoords.size())
            {
              theCoords.push_back(aCoord);
            }
          else
            {
              theCoords[theCoordSize-1] = aCoord;
            }
          theVariable->setValue(theCoordSize);
        }
    }
  void addCompCoord(unsigned int aCoord)
    {
      theLattice[aCoord].id = theID;
      theCompCoords->push_back(aCoord);
      ++theCoordSize;
      theVariable->setValue(theCoordSize);
    }
  unsigned int compCoordSize()
    {
      return theCompCoords->size();
    }
  unsigned int getCompCoord(unsigned int index)
    {
      return (*theCompCoords)[index];
    }
  //it is soft remove because the id of the molecule is not changed:
  void softRemoveCoord(unsigned int aCoord)
    {
      if(!isVacant)
        {
          for(unsigned int i(0); i < theCoordSize; ++i)
            {
              if(theCoords[i] == aCoord)
                {
                  theCoords[i] = theCoords[--theCoordSize];
                  theVariable->setValue(theCoordSize);
                  return;
                }
            }
        }
    }
  void removeCoord(unsigned int aCoord)
    {
      if(!isVacant)
        {
          for(unsigned int i(0); i < theCoordSize; ++i)
            {
              if(theCoords[i] == aCoord)
                {
                  theLattice[aCoord].id = theVacantID;
                  theCoords[i] = theCoords[--theCoordSize];
                  theVariable->setValue(theCoordSize);
                  return;
                }
            }
        }
    }
  //Used to remove all molecules and free memory used to store the molecules
  void clearCoords()
    {
      theCoords.resize(0);
      theCoordSize = 0;
      theVariable->setValue(0);
    }
  //Used by the SpatiocyteStepper when resetting an interation, so must
  //clear the whole compartment using theComp->vacantSpecies->getVacantID():
  void removeCoords()
    {
      if(!isCompVacant)
        {
          for(unsigned int i(0); i < theCoordSize; ++i)
            {
              theLattice[theCoords[i]].id = theVacantSpecies->getID();
            }
          theCoordSize = 0;
          theVariable->setValue(theCoordSize);
        }
    }
  int getPopulateCoordSize()
    {
      return theInitCoordSize-theCoordSize;
    }
  int getInitCoordSize()
    {
      return theInitCoordSize;
    }
  void initMoleculeOrigins()
    {
      theCoordOrigins.resize(theCoordSize);
      for(unsigned int i(0); i < theCoordSize; ++i)
        {
          Origin& anOrigin(theCoordOrigins[i]);
          anOrigin.point = theStepper->coord2point(theCoords[i]);
          anOrigin.row = 0;
          anOrigin.layer = 0;
          anOrigin.col = 0;
        }
    }
  void removeBoundaryCoords()
    {
      for(unsigned int i(0); i < theCoordSize; ++i)
        {
          if(theStepper->isBoundaryCoord(theCoords[i], theDimension))
            {
              std::cout << "is still there" << std::endl;
            }
        }
      theVariable->setValue(theCoordSize);
    }
  void relocateBoundaryCoords()
    {
      for(unsigned int i(0); i < theCoordSize; ++i)
        {
          Origin anOrigin(theCoordOrigins[i]);
          unsigned int periodicCoord(theStepper->getPeriodicCoord(
                                                theCoords[i],
                                                theDimension, &anOrigin));
          if(theLattice[periodicCoord].id == theVacantID)
            {
              theLattice[theCoords[i]].id = theVacantID;
              theCoords[i] = periodicCoord;
              theLattice[theCoords[i]].id = theID;
              theCoordOrigins[i] = anOrigin;
            }
        }
    }
  int getVacantID() const
    {
      return theVacantID;
    }
  Species* getVacantSpecies()
    {
      return theVacantSpecies;
    }
  void setVacantSpecies(Species* aVacantSpecies)
    {
      theVacantSpecies = aVacantSpecies;
      theVacantID = aVacantSpecies->getID();
    }
   const std::vector<double>& getBendAngles() const
    {
      return theBendAngles;
    }
  const Point getWestPoint() const
    {
      Point aWestPoint(theComp->centerPoint);
      aWestPoint.x = theComp->centerPoint.x-theComp->lengthX/2+
        theComp->lengthY/2;
      return aWestPoint;
    }
  const Point getEastPoint() const
    {
      Point anEastPoint(theComp->centerPoint); 
      anEastPoint.x = theComp->centerPoint.x+theComp->lengthX/2-
        theComp->lengthY/2;
      return anEastPoint;
    }
  double getCompRadius() const
    {
      return theComp->lengthY/2;
    }
  double getRadius() const
    {
      return theRadius;
    }
  void setRadius(double aRadius)
    {
      theRadius = aRadius;
    }
  Species* getDiffusionInfluencedReactantPair()
    {
      if(theDiffusionInfluencedReactantPairs.empty())
        {
          return NULL;
        }
      return theDiffusionInfluencedReactantPairs[0];
    }
  double getReactionProbability(int anID)
    {
      return theReactionProbabilities[anID];
    }
  double getMaxReactionProbability()
    {
      double maxProbability(0);
      for(std::vector<double>::const_iterator i(theReactionProbabilities.begin()); 
          i != theReactionProbabilities.end(); ++i)
        {
          if(maxProbability < *i)
            {
              maxProbability = *i;
            }
        }
      return maxProbability;
    }
  void rescaleReactionProbabilities(double aWalkProbability)
    {
      theWalkProbability = aWalkProbability;
      for(std::vector<double>::iterator i(theReactionProbabilities.begin()); 
          i != theReactionProbabilities.end(); ++i)
        {
          *i = (*i)*aWalkProbability;
        }
    }
  void setDiffusionInterval(double anInterval)
    {
      theDiffusionInterval = anInterval;
    }
  unsigned int getRandomCoord()
    {
      if(theCoordSize == 0)
        {
          std::cout << theVariable->getFullID().asString() << std::endl;
          std::cout << "Species size error:" <<
            theVariable->getValue() << std::endl;
        }
      return theCoords[gsl_rng_uniform_int(theRng, theCoordSize)];
    }
  void addInterruptedProcess(SpatiocyteProcessInterface* aProcess)
    {
      theInterruptedProcesses.push_back(aProcess);
    }
  int getBendIndex(double aBendAngle)
    {
      for(unsigned int i(0); i != theBendAngles.size(); ++i)
        {
          if(theBendAngles[i] == aBendAngle)
            {
              return i;
            }
        }
      return 0;
    }
  unsigned int getRandomAdjoiningCoord(unsigned int coord, int searchVacant)
    {
      Voxel& source(theLattice[coord]);
      std::vector<unsigned int> compCoords;
      if(searchVacant)
        { 
          for(unsigned int i(0); i != source.adjoiningSize; ++i)
            {
              unsigned int aCoord(source.adjoiningCoords[i]);
              if(theLattice[aCoord].id == theVacantID)
                {
                  compCoords.push_back(aCoord);
                }
            }
        }
      else
        {
          for(unsigned int i(0); i != source.adjoiningSize; ++i)
            {
              unsigned int aCoord(source.adjoiningCoords[i]);
              if(theStepper->id2Comp(theLattice[aCoord].id) == theComp)
                {
                  compCoords.push_back(aCoord);
                }
            }
        }
      return getRandomVacantCoord(compCoords);
    } 
  unsigned int getBindingSiteAdjoiningCoord(unsigned int coord, int bindingSite)
    {
      Voxel& source(theLattice[coord]);
      if(bindingSite < source.adjoiningSize)
        { 
          unsigned int aCoord(source.adjoiningCoords[bindingSite]);
          if(theLattice[aCoord].id == theVacantID)
            {
              return aCoord;
            }
        }
      return theNullCoord;
    } 
  unsigned int getRandomAdjoiningCoord(unsigned int coord,
                                       Species* aTargetSpecies,
                                       int searchVacant)
    {
      Voxel& source(theLattice[coord]);
      std::vector<unsigned int> compCoords;
      if(searchVacant)
        { 
          for(unsigned int i(0); i != source.adjoiningSize; ++i)
            {
              unsigned int aCoord(source.adjoiningCoords[i]);
              if(theLattice[aCoord].id == aTargetSpecies->getID())
                {
                  compCoords.push_back(aCoord);
                }
            }
        }
      else
        {
          for(unsigned int i(0); i != source.adjoiningSize; ++i)
            {
              unsigned int aCoord(source.adjoiningCoords[i]);
              if(theStepper->id2Comp(theLattice[aCoord].id) == theComp)
                {
                  compCoords.push_back(aCoord);
                }
            }
        }
      return getRandomVacantCoord(compCoords, aTargetSpecies);
    } 
  unsigned int getRandomAdjoiningCoord(unsigned int coord,
                                       unsigned int target,
                                       int searchVacant)
    {
      Voxel& source(theLattice[coord]);
      std::vector<unsigned int> compCoords;
      if(searchVacant)
        { 
          for(unsigned int i(0); i != source.adjoiningSize; ++i)
            {
              unsigned int aCoord(source.adjoiningCoords[i]);
              if(theLattice[aCoord].id == theVacantID && aCoord != target)
                {
                  compCoords.push_back(aCoord);
                }
            }
        }
      else
        {
          for(unsigned int i(0); i != source.adjoiningSize; ++i)
            {
              unsigned int aCoord(source.adjoiningCoords[i]);
              if(theStepper->id2Comp(theLattice[aCoord].id) == theComp &&
                 aCoord != target)
                {
                  compCoords.push_back(aCoord);
                }
            }
        }
      return getRandomVacantCoord(compCoords);
    }
  unsigned int getAdjoiningCoordCnt(unsigned int coord, Species* aTargetSpecies)
    {
      Voxel& source(theLattice[coord]);
      unsigned int cnt(0);
      for(unsigned int i(0); i != source.adjoiningSize; ++i)
        {
          if(theLattice[source.adjoiningCoords[i]].id == 
             aTargetSpecies->getID())
            {
              ++cnt;
            }
        }
      return cnt;
    }
  unsigned int getRandomAdjoiningCoord(unsigned int coord,
                                       unsigned int targetA,
                                       unsigned int targetB, int searchVacant)
    {
      Voxel& source(theLattice[coord]);
      std::vector<unsigned int> compCoords;
      if(searchVacant)
        { 
          for(unsigned int i(0); i != source.adjoiningSize; ++i)
            {
              unsigned int aCoord(source.adjoiningCoords[i]);
              if(theLattice[aCoord].id == theVacantID &&
                 aCoord != targetA && aCoord != targetB)
                {
                  compCoords.push_back(aCoord);
                }
            }
        }
      else
        {
          for(unsigned int i(0); i != source.adjoiningSize; ++i)
            {
              unsigned int aCoord(source.adjoiningCoords[i]);
              if(theStepper->id2Comp(theLattice[aCoord].id) == theComp &&
                 aCoord != targetA && aCoord != targetB)
                {
                  compCoords.push_back(aCoord);
                }
            }
        }
      return getRandomVacantCoord(compCoords);
    }
  unsigned int getRandomVacantCoord(std::vector<unsigned int>& aCoords)
    {
      if(aCoords.size())
        {
          const int r(gsl_rng_uniform_int(theRng, aCoords.size())); 
          unsigned int aCoord(aCoords[r]);
          if(theLattice[aCoord].id == theVacantID)
            {
              return aCoord;
            }
        }
      return theNullCoord;
    }
  unsigned int getRandomVacantCoord(std::vector<unsigned int>& aCoords,
                                    Species* aVacantSpecies)
    {
      if(aCoords.size())
        {
          const int r(gsl_rng_uniform_int(theRng, aCoords.size())); 
          unsigned int aCoord(aCoords[r]);
          if(theLattice[aCoord].id == aVacantSpecies->getID())
            {
              return aCoord;
            }
        }
      return theNullCoord;
    }
  unsigned int getRandomCompCoord(int searchVacant)
    {
      Species* aVacantSpecies(theComp->vacantSpecies);
      int aSize(aVacantSpecies->compCoordSize());
      int r(gsl_rng_uniform_int(theRng, aSize));
      if(searchVacant)
        {
          for(int i(r); i != aSize; ++i)
            {
              unsigned int aCoord(aVacantSpecies->getCompCoord(i));
              if(theLattice[aCoord].id == theVacantID)
                {
                  return aCoord;
                }
            }
          for(int i(0); i != r; ++i)
            {
              unsigned int aCoord(aVacantSpecies->getCompCoord(i));
              if(theLattice[aCoord].id == theVacantID)
                {
                  return aCoord;
                }
            }
        }
      else
        {
          unsigned int aCoord(aVacantSpecies->getCompCoord(r));
          if(theLattice[aCoord].id == theVacantID)
            {
              return aCoord;
            }
        }
      return theNullCoord;
    }
  unsigned int getRandomAdjoiningCompCoord(Comp* aComp, int searchVacant)
    {
      int aSize(theVacantSpecies->size());
      int r(gsl_rng_uniform_int(theRng, aSize)); 
      unsigned int aCoord(theVacantSpecies->getCoord(r));
      return getRandomAdjoiningCoord(aCoord, searchVacant);
    }
  void setLatticeProperties(unsigned int aNullCoord)
    {
      theNullCoord = aNullCoord;
    }
private:
  bool isCentered;
  bool isCompVacant;
  bool isDiffusing;
  bool isDiffusiveVacant;
  bool isFixedAdjoins;
  bool isGaussianPopulation;
  bool isInContact;
  bool isOffLattice;
  bool isPolymer;
  bool isReactiveVacant;
  bool isSubunitInitialized;
  bool isVacant;
  const unsigned short theID;
  unsigned int theCollision;
  unsigned int theDimension;
  unsigned int theInitCoordSize;
  unsigned int theCoordSize;
  unsigned int theNullCoord;
  unsigned int theAdjoiningCoordSize;
  int thePolymerDirectionality;
  int theVacantID;
  double D;
  double theDiffusionInterval;
  double theWalkProbability;
  double theRadius;
  const gsl_rng* theRng;
  Species* theVacantSpecies;
  Comp* theComp;
  MoleculePopulateProcessInterface* thePopulateProcess;
  SpatiocyteStepper* theStepper;
  Variable* theVariable;
  std::vector<bool> theFinalizeReactions;
  std::vector<unsigned int> collisionCnts;
  std::vector<double> theBendAngles;
  std::vector<double> theReactionProbabilities;
  std::vector<unsigned int> theCoords;
  std::vector<unsigned int>* theCompCoords;
  std::vector<Species*> theDiffusionInfluencedReactantPairs;
  std::vector<DiffusionInfluencedReactionProcessInterface*> 
    theDiffusionInfluencedReactions;
  std::vector<SpatiocyteProcessInterface*> theInterruptedProcesses;
  std::vector<Origin> theCoordOrigins;
  std::vector<Voxel>& theLattice;
};


#endif /* __SpatiocyteSpecies_hpp */

