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
#include "SpatiocyteCommon.hpp"
#include "SpatiocyteStepper.hpp"
#include "SpatiocyteProcessInterface.hpp"
#include "DiffusionInfluencedReactionProcessInterface.hpp"
#include "MoleculePopulateProcessInterface.hpp"

// The size of Voxel must be 128 bytes to avoid cacheline splits
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
 * that theMolecules list is not updated when other species diffuse on it. 
 * There are four possible types of vacant species:
 * 1. !isDiffusiveVacant && !isReactiveVacant: a vacant species
 *    whose theMolecules list and theMoleculeSize are never updated. This is
 *    the most basic version. This vacant species never diffuses, or reacts
 *    using SNRP. 
 * 2. !isDiffusiveVacant && isReactiveVacant: a vacant species which is a
 *    substrate of a SNRP reaction. In this case theMoleculeSize is updated
 *    before the step interval of SNRP is calculated, and theMolecules list is
 *    updated before the SNRP reaction (fire) is executed to get a valid
 *    molecule list. Set by SNRP.
 * 3. isDiffusiveVacant && !isReactiveVacant: a vacant species which also
 *    diffuses. In this case, theMolecules list and theMoleculesSize are
 *    updated just before it is diffused. Set by DiffusionProcess.   
 * 4. isDiffusiveVacant && isReactiveVacant: a vacant species which reacts
 *    using SNRP and also diffuses.
 * 5. isCompVacant: the VACANT species declared for each compartment. It can
 *    also be isDiffusiveVacant and isReactiveVacant. Set during compartment
 *    registration. It also persistently stores all the compartment voxels.
 *    Referred to as theVacantSpecies. For isCompVacant, initially there
 *    are no molecules in its list. All voxels are stored in theCompVoxels. Only
 *    if it is updated before being called by VisualizationLogProcess, SNRP or
 *    DiffusionProcess, it will be theMolecules will be populated with the
 *    CompVoxels.
 * 6. isVacant {isCompVacant; isDiffusiveVacant; isReactiveVacant): the
 *    general name used to identify either isCompVacant, isDiffusiveVacant or
 *    isReactiveVacant to reduce comparison operations.
 */

class Species
{
public:
  Species(SpatiocyteStepper* aStepper, Variable* aVariable, int anID, 
          int anInitMoleculeSize, const gsl_rng* aRng, double voxelRadius):
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
    theInitMoleculeSize(anInitMoleculeSize),
    theMoleculeSize(0),
    D(0),
    theDiffusionInterval(libecs::INF),
    theWalkProbability(1),
    theRadius(voxelRadius),
    theRng(aRng),
    thePopulateProcess(NULL),
    theStepper(aStepper),
    theVariable(aVariable),
    theCompVoxels(&theMolecules) {}
  ~Species() {}
  void initialize(int speciesSize, int anAdjoiningVoxelSize)
    {
      theAdjoiningVoxelSize = anAdjoiningVoxelSize;
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
      else if(theMoleculeSize)
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
      else if(theMoleculeSize)
        {
          std::cout << "Species:" << theVariable->getFullID().asString() <<
            " not MoleculePopulated." << std::endl;
        }
    }
  void populateCompUniformSparse()
    {
      if(thePopulateProcess)
        {
          thePopulateProcess->populateUniformSparse(this);
        }
      else if(theMoleculeSize)
        {
          std::cout << "Species:" << theVariable->getFullID().asString() <<
            " not MoleculePopulated." << std::endl;
        }
    }
  void populateUniformOnDiffusiveVacant()
    {
      if(thePopulateProcess)
        {
          thePopulateProcess->populateUniformOnDiffusiveVacant(this);
        }
      else if(theMoleculeSize)
        {
          std::cout << "Species:" << theVariable->getFullID().asString() <<
            " not MoleculePopulated." << std::endl;
        }
    }
  Variable* getVariable() const
    {
      return theVariable;
    }
  std::vector<unsigned int> getSourceCoords()
    {
      std::vector<unsigned int> aCoords;
      for(unsigned int i(0); i != theMoleculeSize; ++i)
        {
          std::vector<Voxel*>& 
            aSourceVoxels(theMolecules[i]->subunit->sourceVoxels);
          for(unsigned int j(0); j != aSourceVoxels.size(); ++j)
            {
              if(aSourceVoxels[j])
                {
                  aCoords.push_back(aSourceVoxels[j]->coord);
                }
            }
        }
      return aCoords;
    }
  std::vector<unsigned int> getTargetCoords()
    {
      std::vector<unsigned int> aCoords;
      for(unsigned int i(0); i != theMoleculeSize; ++i)
        {
          std::vector<Voxel*>& 
            aTargetVoxels(theMolecules[i]->subunit->targetVoxels);
          for(unsigned int j(0); j != aTargetVoxels.size(); ++j)
            {
              if(aTargetVoxels[j])
                {
                  aCoords.push_back(aTargetVoxels[j]->coord);
                }
            }
        }
      return aCoords;
    }
  std::vector<unsigned int> getSharedCoords()
    {
      std::vector<unsigned int> aCoords;
      for(unsigned int i(0); i != theMoleculeSize; ++i)
        {
          std::vector<Voxel*>& 
            aSharedLipids(theMolecules[i]->subunit->sharedLipids);
          for(unsigned int j(0); j != aSharedLipids.size(); ++j)
            {
              if(aSharedLipids[j])
                {
                  aCoords.push_back(aSharedLipids[j]->coord);
                }
            }
        }
      return aCoords;
    }
  std::vector<Voxel*>& getMolecules()
    {
      return theMolecules;
    }
  unsigned int size() const
    {
      return theMoleculeSize;
    }
  unsigned int getCoord(int anIndex)
    {
      return theMolecules[anIndex]->coord;
    }
  Point getPoint(int anIndex)
    {
      if(isOffLattice)
        {
          return *theMolecules[anIndex]->point;
        }
      else if(isPolymer)
        {
          return theMolecules[anIndex]->subunit->subunitPoint;
        }
      return theStepper->coord2point(theMolecules[anIndex]->coord);
    }
  Voxel* getMolecule(int anIndex)
    {
      return theMolecules[anIndex];
    }
  unsigned short getID() const
    {
      return theID;
    }
  double getMeanSquaredDisplacement()
    {
      if(!theMoleculeSize)
        {
          return 0;
        }
      double aDisplacement(0);
      for(unsigned int i(0); i < theMoleculeSize; ++i)
        {
          Point aCurrentPoint(theStepper->getPeriodicPoint(
                                                 theMolecules[i]->coord,
                                                 theDimension,
                                                 &theMoleculeOrigins[i]));
          double aDistance(getDistance(&theMoleculeOrigins[i].point,
                                       &aCurrentPoint));
          aDisplacement += aDistance*aDistance;
        }
      return
        aDisplacement*pow(theRadius*2, 2)/theMoleculeSize;
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
      theInitMoleculeSize = theMoleculeSize;
      getVariable()->setValue(theMoleculeSize);
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
      return isCompVacant;
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
      return theMoleculeSize == theInitMoleculeSize;
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
  void walk()
    {
      for(unsigned int i(0); i < theMoleculeSize; ++i)
        {
          Voxel* source(theMolecules[i]);
          int size;
          if(isFixedAdjoins)
            {
              size = theAdjoiningVoxelSize;
            }
          else
            {
              size = source->diffuseSize;
            }
          Voxel* target(source->adjoiningVoxels[
                        gsl_rng_uniform_int(theRng, size)]);
          /*
          if(source == target)
            {
              std::cout << "SpatiocyteSpecies source=target error" << std::endl;
            }
            */
          if(target->id == theVacantID)
            {
              if(theWalkProbability == 1 ||
                 gsl_rng_uniform(theRng) < theWalkProbability)
                {
                  target->id = theID;
                  source->id = theVacantID;
                  theMolecules[i] = target;
                }
            }
          else if(theDiffusionInfluencedReactions[target->id] != NULL)
            {
              //If it meets the reaction probability:
              if(gsl_rng_uniform(theRng) < theReactionProbabilities[target->id])
                { 
                  Species* targetSpecies(theStepper->id2species(target->id));
                  DiffusionInfluencedReactionProcessInterface* aReaction(
                             theDiffusionInfluencedReactions[target->id]);
                  if(aReaction->react(source, target))
                    {
                      //Soft remove the source molecule, i.e.,
                      //keep the id intact:
                      theMolecules[i--] = theMolecules[--theMoleculeSize];
                      theVariable->setValue(theMoleculeSize);
                      //Soft remove the target molecule:
                      targetSpecies->softRemoveMolecule(target);
                      theFinalizeReactions[targetSpecies->getID()] = true;
                    }
                }
            }
        }
    }
  void walkVacant()
    {
      updateVacantMolecules();
      for(unsigned int i(0); i < theMoleculeSize; ++i)
        {
          Voxel* source(theMolecules[i]);
          int size;
          if(isFixedAdjoins)
            {
              size = theAdjoiningVoxelSize;
            }
          else
            {
              size = source->diffuseSize;
            }
          Voxel* target(source->adjoiningVoxels[
                        gsl_rng_uniform_int(theRng, size)]);
          if(target->id == theVacantID)
            {
              if(theWalkProbability == 1 ||
                 gsl_rng_uniform(theRng) < theWalkProbability)
                {
                  target->id = theID;
                  source->id = theVacantID;
                }
            }
          /*
          else if(theDiffusionInfluencedReactions[target->id] != NULL)
            {
              //If it meets the reaction probability:
              if(gsl_rng_uniform(theRng) <
                 theReactionProbabilities[target->id])
                { 
                  Species* targetSpecies(theStepper->id2species(target->id));
                  DiffusionInfluencedReactionProcessInterface* aReaction(
                             theDiffusionInfluencedReactions[target->id]);
                  if(aReaction->react(source, target))
                    {
                      //Soft remove the target molecule:
                      targetSpecies->softRemoveMolecule(target);
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
      int newMoleculeSize(0);
      for(unsigned int i(0); i < theMoleculeSize; ++i) 
        {
          Voxel* aVoxel(theMolecules[i]);
          if(theStepper->isRemovableEdgeCoord(aVoxel->coord, theComp))
            {
              Comp* aSuperComp(
                 theStepper->system2Comp(theComp->system->getSuperSystem())); 
              aSuperComp->vacantSpecies->addCompVoxel(aVoxel);
            }
          else 
            { 
              theMolecules[newMoleculeSize] = aVoxel; 
              ++newMoleculeSize; 
            }
        }
      theMoleculeSize = newMoleculeSize;
      theVariable->setValue(theMoleculeSize);
    }
  void removePeriodicEdgeVoxels()
    {
      int newMoleculeSize(0);
      for(unsigned int i(0); i < theMoleculeSize; ++i) 
        {
          Voxel* aVoxel(theMolecules[i]);
          if(theStepper->isPeriodicEdgeCoord(aVoxel->coord, theComp))
            {
              aVoxel->id = theStepper->getNullID();
            }
          else 
            { 
              theMolecules[newMoleculeSize] = aVoxel; 
              ++newMoleculeSize; 
            }
        }
      theMoleculeSize = newMoleculeSize;
      theVariable->setValue(theMoleculeSize);
    }
  void updateSpecies()
    {
      if(isCompVacant && (isDiffusiveVacant || isReactiveVacant))
        {
          theCompVoxels = new std::vector<Voxel*>;
          for(unsigned int i(0); i != theMoleculeSize; ++i)
            { 
              theCompVoxels->push_back(theMolecules[i]);
            }
        }
    }
  //If it isReactiveVacant it will only be called by SNRP when it is substrate
  //If it isDiffusiveVacant it will only be called by DiffusionProcess before
  //being diffused. So we need to only check if it isVacant:
  void updateMolecules()
    {
      if(isDiffusiveVacant || isReactiveVacant)
        {
          updateVacantMolecules();
        }
    }
  //If it isReactiveVacant it will only be called by SNRP when it is substrate:
  void updateMoleculeSize()
    {
      if(isDiffusiveVacant || isReactiveVacant)
        {
          updateVacantMoleculeSize();
        }
    }
  //Even if it is a isCompVacant, this method will be called by
  //VisualizationLogProcess, or SNRP if it is Reactive or DiffusionProcess
  //if it is Diffusive:
  void updateVacantMolecules()
    {
      theMoleculeSize = 0;
      int aSize(theVacantSpecies->compVoxelSize());
      for(int i(0); i != aSize; ++i)
        { 
          Voxel* aMolecule(theVacantSpecies->getCompVoxel(i));
          if(aMolecule->id == theID)
            {
              ++theMoleculeSize;
              if(theMoleculeSize > theMolecules.size())
                {
                  theMolecules.push_back(aMolecule);
                }
              else
                {
                  theMolecules[theMoleculeSize-1] = aMolecule;
                }
            }
        }
      theVariable->setValue(theMoleculeSize);
    }
  void updateVacantMoleculeSize()
    {
      theMoleculeSize = 0;
      int aSize(theVacantSpecies->compVoxelSize());
      for(int i(0); i != aSize; ++i)
        { 
          Voxel* aMolecule(theVacantSpecies->getCompVoxel(i));
          if(aMolecule->id == theID)
            {
              ++theMoleculeSize;
            }
        }
      if(theMoleculeSize > theMolecules.size())
        {
          theMolecules.resize(theMoleculeSize);
        }
      theVariable->setValue(theMoleculeSize);
    }
  void addMolecule(Voxel* aMolecule)
    {
      aMolecule->id = theID;
      if(!isVacant)
        {
          ++theMoleculeSize;
          if(theMoleculeSize > theMolecules.size())
            {
              theMolecules.push_back(aMolecule);
            }
          else
            {
              theMolecules[theMoleculeSize-1] = aMolecule;
            }
          theVariable->setValue(theMoleculeSize);
        }
    }
  void addCompVoxel(Voxel* aVoxel)
    {
      aVoxel->id = theID;
      theCompVoxels->push_back(aVoxel);
      ++theMoleculeSize;
      theVariable->setValue(theMoleculeSize);
    }
  unsigned int compVoxelSize()
    {
      return theCompVoxels->size();
    }
  Voxel* getCompVoxel(unsigned int index)
    {
      return (*theCompVoxels)[index];
    }
  //it is soft remove because the id of the molecule is not changed:
  void softRemoveMolecule(Voxel* aMolecule)
    {
      if(!isVacant)
        {
          for(unsigned int i(0); i < theMoleculeSize; ++i)
            {
              if(theMolecules[i] == aMolecule)
                {
                  theMolecules[i] = theMolecules[--theMoleculeSize];
                  theVariable->setValue(theMoleculeSize);
                  return;
                }
            }
        }
    }
  void removeMolecule(Voxel* aMolecule)
    {
      if(!isVacant)
        {
          for(unsigned int i(0); i < theMoleculeSize; ++i)
            {
              if(theMolecules[i] == aMolecule)
                {
                  aMolecule->id = theVacantID;
                  theMolecules[i] = theMolecules[--theMoleculeSize];
                  theVariable->setValue(theMoleculeSize);
                  return;
                }
            }
        }
    }
  //Used to remove all molecules and free memory used to store the molecules
  void clearMolecules()
    {
      theMolecules.resize(0);
      theMoleculeSize = 0;
      theVariable->setValue(0);
    }
  //Used by the SpatiocyteStepper when resetting an interation, so must
  //clear the whole compartment using theComp->vacantSpecies->getVacantID():
  void removeMolecules()
    {
      for(unsigned int i(0); i < theMoleculeSize; ++i)
        {
          theMolecules[i]->id = theVacantSpecies->getID();
        }
      theMoleculeSize = 0;
      theVariable->setValue(theMoleculeSize);
    }
  int getPopulateMoleculeSize()
    {
      return theInitMoleculeSize-theMoleculeSize;
    }
  int getInitMoleculeSize()
    {
      return theInitMoleculeSize;
    }
  void initMoleculeOrigins()
    {
      theMoleculeOrigins.resize(theMoleculeSize);
      for(unsigned int i(0); i < theMoleculeSize; ++i)
        {
          Origin& anOrigin(theMoleculeOrigins[i]);
          anOrigin.point = theStepper->coord2point(theMolecules[i]->coord);
          anOrigin.row = 0;
          anOrigin.layer = 0;
          anOrigin.col = 0;
        }
    }
  void removeBoundaryMolecules()
    {
      for(unsigned int i(0); i < theMoleculeSize; ++i)
        {
          if(theStepper->isBoundaryCoord(
                               theMolecules[i]->coord, theDimension))
            {
              std::cout << "is still there" << std::endl;
            }
        }
      theVariable->setValue(theMoleculeSize);
    }
  void relocateBoundaryMolecules()
    {
      for(unsigned int i(0); i < theMoleculeSize; ++i)
        {
          Origin anOrigin(theMoleculeOrigins[i]);
          Voxel* periodicVoxel(theStepper->getPeriodicVoxel(
                                                theMolecules[i]->coord,
                                                theDimension, &anOrigin));
          if(periodicVoxel != NULL && 
             periodicVoxel->id == theVacantID)
            {
              theMolecules[i]->id = theVacantID;
              theMolecules[i] = periodicVoxel;
              theMolecules[i]->id = theID;
              theMoleculeOrigins[i] = anOrigin;
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
  Voxel* getRandomMolecule()
    {
      if(theMoleculeSize == 0)
        {
          std::cout << theVariable->getFullID().asString() << std::endl;
          std::cout << "Species size error:" <<
            theVariable->getValue() << std::endl;
        }
      return theMolecules[gsl_rng_uniform_int(theRng, theMoleculeSize)];
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
  Voxel* getRandomAdjoiningVoxel(Voxel* source)
    {
      std::vector<Voxel*> CompVoxels;
      if(theStepper->getSearchVacant())
        { 
          for(unsigned int i(0); i != source->adjoiningSize; ++i)
            {
              Voxel* aVoxel(source->adjoiningVoxels[i]);
              if(aVoxel->id == theVacantID)
                {
                  CompVoxels.push_back(aVoxel);
                }
            }
        }
      else
        {
          for(unsigned int i(0); i != source->adjoiningSize; ++i)
            {
              Voxel* aVoxel(source->adjoiningVoxels[i]);
              if(theStepper->id2Comp(aVoxel->id) == theComp)
                {
                  CompVoxels.push_back(aVoxel);
                }
            }
        }
      return getRandomVacantVoxel(&CompVoxels);
    } 
  Voxel* getRandomAdjoiningVoxel(Voxel* source, int bindingSite)
    {
      if(bindingSite < source->adjoiningSize)
        { 
          Voxel* aVoxel(source->adjoiningVoxels[bindingSite]);
          if(aVoxel->id == theVacantID)
            {
              return aVoxel;
            }
        }
      return NULL;
    } 
  Voxel* getRandomAdjoiningVoxel(Voxel* source, Species* aVacantSpecies)
    {
      std::vector<Voxel*> CompVoxels;
      if(theStepper->getSearchVacant())
        { 
          for(unsigned int i(0); i != source->adjoiningSize; ++i)
            {
              Voxel* aVoxel(source->adjoiningVoxels[i]);
              if(aVoxel->id == aVacantSpecies->getID())
                {
                  CompVoxels.push_back(aVoxel);
                }
            }
        }
      else
        {
          for(unsigned int i(0); i != source->adjoiningSize; ++i)
            {
              Voxel* aVoxel(source->adjoiningVoxels[i]);
              if(theStepper->id2Comp(aVoxel->id) == theComp)
                {
                  CompVoxels.push_back(aVoxel);
                }
            }
        }
      return getRandomVacantVoxel(&CompVoxels, aVacantSpecies);
    } 
  Voxel* getRandomAdjoiningVoxel(Voxel* source, Voxel* target)
    {
      std::vector<Voxel*> CompVoxels;
      if(theStepper->getSearchVacant())
        { 
          for(unsigned int i(0); i != source->adjoiningSize; ++i)
            {
              Voxel* aVoxel(source->adjoiningVoxels[i]);
              if(aVoxel->id == theVacantID && aVoxel != target)
                {
                  CompVoxels.push_back(aVoxel);
                }
            }
        }
      else
        {
          for(unsigned int i(0); i != source->adjoiningSize; ++i)
            {
              Voxel* aVoxel(source->adjoiningVoxels[i]);
              if(theStepper->id2Comp(aVoxel->id) == theComp &&
                 aVoxel != target)
                {
                  CompVoxels.push_back(aVoxel);
                }
            }
        }
      return getRandomVacantVoxel(&CompVoxels);
    }
  Voxel* getRandomAdjoiningVoxel(Voxel* source, Voxel* targetA, Voxel* targetB)
    {
      std::vector<Voxel*> CompVoxels;
      if(theStepper->getSearchVacant())
        { 
          for(unsigned int i(0); i != source->adjoiningSize; ++i)
            {
              Voxel* aVoxel(source->adjoiningVoxels[i]);
              if(aVoxel->id == theVacantID &&
                 aVoxel != targetA && aVoxel != targetB)
                {
                  CompVoxels.push_back(aVoxel);
                }
            }
        }
      else
        {
          for(unsigned int i(0); i != source->adjoiningSize; ++i)
            {
              Voxel* aVoxel(source->adjoiningVoxels[i]);
              if(theStepper->id2Comp(aVoxel->id) == theComp &&
                 aVoxel != targetA && aVoxel != targetB)
                {
                  CompVoxels.push_back(aVoxel);
                }
            }
        }
      return getRandomVacantVoxel(&CompVoxels);
    }
  Voxel* getRandomVacantVoxel(std::vector<Voxel*>* voxels)
    {
      if(voxels->size())
        {
          const int r(gsl_rng_uniform_int(theRng, voxels->size())); 
          Voxel* aVoxel((*voxels)[r]);
          if(aVoxel->id == theVacantID)
            {
              return aVoxel;
            }
        }
      return NULL;
    }
  Voxel* getRandomVacantVoxel(std::vector<Voxel*>* voxels,
                              Species* aVacantSpecies)
    {
      if(voxels->size())
        {
          const int r(gsl_rng_uniform_int(theRng, voxels->size())); 
          Voxel* aVoxel((*voxels)[r]);
          if(aVoxel->id == aVacantSpecies->getID())
            {
              return aVoxel;
            }
        }
      return NULL;
    }
  Voxel* getRandomCompVoxel()
    {
      Species* aVacantSpecies(theComp->vacantSpecies);
      int aSize(aVacantSpecies->compVoxelSize());
      int r(gsl_rng_uniform_int(theRng, aSize));
      if(theStepper->getSearchVacant())
        {
          for(int i(r); i != aSize; ++i)
            {
              Voxel* aVoxel(aVacantSpecies->getCompVoxel(i));
              if(aVoxel->id == theVacantID)
                {
                  return aVoxel;
                }
            }
          for(int i(0); i != r; ++i)
            {
              Voxel* aVoxel(aVacantSpecies->getCompVoxel(i));
              if(aVoxel->id == theVacantID)
                {
                  return aVoxel;
                }
            }
        }
      else
        {
          Voxel* aVoxel(aVacantSpecies->getCompVoxel(r));
          if(aVoxel->id == theVacantID)
            {
              return aVoxel;
            }
        }
      return NULL;
    }
  Voxel* getRandomAdjoiningCompVoxel(Comp* aComp)
    {
      int aSize(theVacantSpecies->size());
      int r(gsl_rng_uniform_int(theRng, aSize)); 
      Voxel* aVoxel(theVacantSpecies->getMolecule(r));
      return getRandomAdjoiningVoxel(aVoxel);
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
  unsigned int theDimension;
  unsigned int theInitMoleculeSize;
  unsigned int theMoleculeSize;
  unsigned int theAdjoiningVoxelSize;
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
  std::vector<double> theBendAngles;
  std::vector<double> theReactionProbabilities;
  std::vector<Voxel*> theMolecules;
  std::vector<Voxel*>* theCompVoxels;
  std::vector<Species*> theDiffusionInfluencedReactantPairs;
  std::vector<DiffusionInfluencedReactionProcessInterface*> 
    theDiffusionInfluencedReactions;
  std::vector<SpatiocyteProcessInterface*> theInterruptedProcesses;
  std::vector<Origin> theMoleculeOrigins;
};


#endif /* __SpatiocyteSpecies_hpp */

