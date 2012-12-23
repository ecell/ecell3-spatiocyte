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
#include "SpatiocyteNextReactionProcess.hpp"
#include "DiffusionInfluencedReactionProcess.hpp"
#include "MoleculePopulateProcessInterface.hpp"

static double getDistance(Point& aSourcePoint, Point& aDestPoint)
{
  return sqrt(pow(aDestPoint.x-aSourcePoint.x, 2)+
              pow(aDestPoint.y-aSourcePoint.y, 2)+
              pow(aDestPoint.z-aSourcePoint.z, 2));
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
 *    whose theCoords list and theMoleculeSize are never updated. This is
 *    the most basic version. This vacant species never diffuses, or reacts
 *    using SNRP. 
 * 2. !isDiffusiveVacant && isReactiveVacant: a vacant species which is a
 *    substrate of a SNRP reaction. In this case theMoleculeSize is updated
 *    before the step interval of SNRP is calculated, and theCoords list is
 *    updated before the SNRP reaction (fire) is executed to get a valid
 *    molecule list. Set by SNRP.
 * 3. isDiffusiveVacant && !isReactiveVacant: a vacant species which also
 *    diffuses. In this case, theCoords list and theMoleculeSize are
 *    updated just before it is diffused. Set by DiffusionProcess.   
 * 4. isDiffusiveVacant && isReactiveVacant: a vacant species which reacts
 *    using SNRP and also diffuses.
 * 5. isCompVacant: the VACANT species declared for each compartment. It can
 *    also be isDiffusiveVacant and isReactiveVacant. Set during compartment
 *    registration. It also persistently stores all the compartment voxels.
 *    Referred to as theVacantSpecies. For isCompVacant, initially there
 *    are no molecules in its list. All voxels are stored in theCompVoxels. Only
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
  Species(SpatiocyteStepper* aStepper, Variable* aVariable,
          unsigned short anID, int anInitCoordSize, const gsl_rng* aRng,
          double voxelRadius, std::vector<Voxel>& aLattice):
    isCentered(false),
    isCompVacant(false),
    isDiffusing(false),
    isDiffusiveVacant(false),
    isFixedAdjoins(false),
    isGaussianPopulation(false),
    isInContact(false),
    isMultiscale(false),
    isOffLattice(false),
    isPolymer(false),
    isReactiveVacant(false),
    isSubunitInitialized(false),
    isTag(false),
    isTagged(false),
    isVacant(false),
    theID(anID),
    theCollision(0),
    theInitCoordSize(anInitCoordSize),
    theMoleculeSize(0),
    D(0),
    theDiffuseRadius(voxelRadius),
    theDiffusionInterval(libecs::INF),
    theMoleculeRadius(voxelRadius),
    theVoxelRadius(voxelRadius),
    theWalkProbability(1),
    theRng(aRng),
    thePopulateProcess(NULL),
    theStepper(aStepper),
    theVariable(aVariable),
    theCompCoords(&theCoords),
    theLattice(aLattice) {}
  ~Species() {}
  void initialize(int speciesSize, int anAdjoinSize,
                  unsigned aNullCoord, unsigned aNullID)
    {
      theAdjoinSize = anAdjoinSize;
      theNullCoord = aNullCoord;
      theNullID = aNullID;
      theSpeciesSize = speciesSize;
      theReactionProbabilities.resize(speciesSize);
      theDiffusionInfluencedReactions.resize(speciesSize);
      theFinalizeReactions.resize(speciesSize);
      theMultiscaleBindIDs.resize(speciesSize);
      theMultiscaleUnbindIDs.resize(speciesSize);
      for(int i(0); i != speciesSize; ++i)
        {
          theDiffusionInfluencedReactions[i] = NULL;
          theReactionProbabilities[i] = 0;
          theFinalizeReactions[i] = false;
        }
      if(theComp)
        {
          setVacantSpecies(theComp->vacantSpecies);
        }
      theNullTag.origin = theNullCoord;
      theNullTag.id = theNullID;
    }
  void setDiffusionInfluencedReaction(
                                    DiffusionInfluencedReactionProcess*
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
          std::cout << "Warning: Species " << getIDString() <<
            " not populated." << std::endl;
        }
    }
  void addTaggedSpecies(Species* aSpecies)
    {
      isTag = true;
      //If one of the tagged species is off-lattice then
      //make the tag species off-lattice. The getPoint method
      //will convert the coord of lattice species to point
      //when logging:
      if(aSpecies->getIsOffLattice())
        {
          isOffLattice = true;
        }
      theTaggedSpeciesList.push_back(aSpecies);
    }
  void addTagSpecies(Species* aSpecies)
    {
      isTagged = true;
      theTagSpeciesList.push_back(aSpecies);
      aSpecies->addTaggedSpecies(this);
    }
  bool getIsTagged()
    {
      return isTagged;
    }
  bool getIsTag()
    {
      return isTag;
    }
  bool getIsPopulateSpecies()
    {
      return (thePopulateProcess != NULL);
    }
  void populateCompUniform(unsigned* voxelIDs, unsigned* aCount)
    {
      if(thePopulateProcess)
        {
          thePopulateProcess->populateUniformDense(this, voxelIDs, aCount);
        }
      else if(theMoleculeSize)
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
      else if(theMoleculeSize)
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
      else if(theMoleculeSize)
        {
          std::cout << "Species:" << theVariable->getFullID().asString() <<
            " not CoordPopulated." << std::endl;
        }
    }
  void populateUniformOnMultiscale()
    {
      if(thePopulateProcess)
        {
          thePopulateProcess->populateUniformOnMultiscale(this);
        }
      else if(theMoleculeSize)
        {
          std::cout << "Species:" << theVariable->getFullID().asString() <<
            " not CoordPopulated." << std::endl;
        }
    }
  Variable* getVariable() const
    {
      return theVariable;
    }
  unsigned size() const
    {
      return theMoleculeSize;
    }
  unsigned getCoord(unsigned anIndex)
    {
      return theCoords[anIndex];
    }
  Point getPoint(unsigned anIndex)
    {
      if(isOffLattice)
        {
          if(theLattice[theCoords[anIndex]].point)
            {
              return *theLattice[theCoords[anIndex]].point;
            }
          return theStepper->coord2point(getCoord(anIndex));
        }
      /*
      else if(isPolymer)
        {
          return theMolecules[anIndex]->subunit->subunitPoint;
        }
        */
      return theStepper->coord2point(getCoord(anIndex));
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
      for(unsigned i(0); i < theMoleculeSize; ++i)
        {
          Point aCurrentPoint(theStepper->getPeriodicPoint(
                                                 getCoord(i),
                                                 theDimension,
                                                 &theMoleculeOrigins[i]));
          double aDistance(getDistance(theMoleculeOrigins[i].point,
                                       aCurrentPoint));
          aDisplacement += aDistance*aDistance;
        }
      return aDisplacement*pow(theDiffuseRadius*2, 2)/theMoleculeSize;
    }
  void setCollision(unsigned aCollision)
    {
      theCollision = aCollision;
    }
  void setIsSubunitInitialized()
    {
      isSubunitInitialized = true;
    }
  void setIsMultiscale()
    {
      isMultiscale = true;
    }
  bool getIsMultiscale()
    {
      return isMultiscale;
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
      theInitCoordSize = theMoleculeSize;
      getVariable()->setValue(theMoleculeSize);
    }
  void finalizeSpecies()
    {
      if(theCollision)
        {
          collisionCnts.resize(theMoleculeSize);
          for(std::vector<unsigned>::iterator 
              i(collisionCnts.begin()); i != collisionCnts.end(); ++i)
            {
              *i = 0;
            }
        }
      //need to shuffle molecules of the compVacant species if it has
      //diffusing vacant species to avoid bias when random walking:
      if(isCompVacant)
        {
          for(unsigned i(0); i != theComp->species.size(); ++i)
            {
              if(theComp->species[i]->getIsDiffusiveVacant())
                {
                  std::random_shuffle(theCoords.begin(), theCoords.end());
                  break;
                }
            }
          if(isTagged)
            {
              for(unsigned i(0); i != theMoleculeSize; ++i)
                {
                  theTags[i].origin = getCoord(i);
                }
            }
        }
    }
  unsigned getCollisionCnt(unsigned anIndex)
    {
      return collisionCnts[anIndex];
    }
  unsigned getCollision() const
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
    }
  void resetFixedAdjoins()
    {
      isFixedAdjoins = false;
      for(unsigned i(0); i != theComp->species.size(); ++i)
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
      return theMoleculeSize == theInitCoordSize;
    }
  double getDiffusionInterval() const
    {
      return theDiffusionInterval;
    }
  unsigned getDimension()
    {
      return theDimension;
    }
  void setDimension(unsigned aDimension)
    {
      theDimension = aDimension;
      if(theDimension == 3)
        {
          isFixedAdjoins = true;
        }
    }
  void resetFinalizeReactions()
    {
      for(unsigned i(0); i != theFinalizeReactions.size(); ++i)
        {
          theFinalizeReactions[i] = false;
        }
    }
  void finalizeReactions()
    {
      for(unsigned i(0); i != theFinalizeReactions.size(); ++i)
        {
          if(theFinalizeReactions[i])
            {
              theDiffusionInfluencedReactions[i]->finalizeReaction();
            }
        }
      for(unsigned i(0); i != theInterruptedProcesses.size(); ++i)
        {
          theInterruptedProcesses[i
            ]->substrateValueChanged(theStepper->getCurrentTime());
        }
    }
  void addCollision(unsigned aCoord)
    {
      for(unsigned i(0); i < theMoleculeSize; ++i)
        {
          if(theCoords[i] == aCoord)
            {
              ++collisionCnts[i];
              return;
            }
        }
      std::cout << "error in species add collision" << std::endl;
    }
  void walk()
    {
      const unsigned beginMoleculeSize(theMoleculeSize);
      unsigned size(theAdjoinSize);
      for(unsigned i(0); i < beginMoleculeSize && i < theMoleculeSize; ++i)
        {
          Voxel& source(theLattice[theCoords[i]]);
          if(!isFixedAdjoins)
            {
              size = source.diffuseSize;
            }
          unsigned tarCoord(source.adjoins[gsl_rng_uniform_int(theRng, size)]);
          unsigned short& targetID(theLattice[tarCoord].id);
          if(targetID == theVacantID)
            {
              if(theWalkProbability == 1 ||
                 gsl_rng_uniform(theRng) < theWalkProbability)
                {
                  targetID = theID;
                  source.id = theVacantID;
                  theCoords[i] = tarCoord;
                }
            }
          else
            {
              if(targetID == theComp->interfaceID)
                {
                  unsigned diffuseSize(theLattice[tarCoord].diffuseSize);
                  unsigned range(theLattice[tarCoord].adjoinSize-diffuseSize);
                  unsigned index(gsl_rng_uniform_int(theRng, range));
                  tarCoord = theLattice[tarCoord].adjoins[diffuseSize+index];
                }
              unsigned short& targetID(theLattice[tarCoord].id);
              if(theDiffusionInfluencedReactions[targetID])
                {
                  //If it meets the reaction probability:
                  if(gsl_rng_uniform(theRng) < 
                     theReactionProbabilities[targetID])
                    { 
                      Species* targetSpecies(theStepper->id2species(targetID));
                      unsigned targetIndex(0);
                      if(targetSpecies->getIsMultiscale() && theVacantSpecies ==
                         targetSpecies->getMultiscaleVacantSpecies())
                        {
                          //Set an invalid index if the target molecule is
                          //an implicitly represented multiscale molecule:
                          targetIndex = targetSpecies->size();
                        }
                      else
                        { 
                          targetIndex = targetSpecies->getMoleculeIndex(tarCoord);
                        }
                      if(theCollision)
                        { 
                          ++collisionCnts[i];
                          targetSpecies->addCollision(tarCoord);
                          if(theCollision != 2)
                            {
                              return;
                            }
                        }
                      unsigned aMoleculeSize(theMoleculeSize);
                      react(theCoords[i], tarCoord, i, targetIndex, targetSpecies);
                      //If the reaction is successful, the last molecule of this
                      //species will replace the pointer of i, so we need to 
                      //decrement i to perform the diffusion on it. However, if
                      //theMoleculeSize didn't decrease, that means the
                      //currently walked molecule was a product of this
                      //reaction and so we don't need to walk it again by
                      //decrementing i.
                      if(theMoleculeSize < aMoleculeSize)
                        {
                          --i;
                        }
                    }
                }
            }
        }
    }
  void walkMultiscale()
    {
      unsigned beginMoleculeSize(theMoleculeSize);
      for(unsigned i(0); i < beginMoleculeSize && i < theMoleculeSize; ++i)
        {
          unsigned srcCoord(theCoords[i]);
          Voxel& source(theLattice[srcCoord]);
          int size(source.diffuseSize);
          unsigned tarCoord(source.adjoins[gsl_rng_uniform_int(theRng, size)]);
          Voxel& target(theLattice[tarCoord]);
          if(target.id == theVacantID)
            {
              if(theWalkProbability == 1 ||
                 gsl_rng_uniform(theRng) < theWalkProbability)
                {
                  if(!isIntersectMultiscale(srcCoord, tarCoord))
                    {
                      removeMultiscaleMolecule(srcCoord);
                      addMultiscaleMolecule(tarCoord);
                      target.id = theID;
                      source.id = theVacantID;
                      theCoords[i] = tarCoord;
                    }
                }
            }
        }
    }
  void walkVacant()
    {
      updateVacantMolecules();
      for(unsigned i(0); i < theMoleculeSize; ++i)
        {
          Voxel* source(&theLattice[theCoords[i]]);
          int size;
          if(isFixedAdjoins)
            {
              size = theAdjoinSize;
            }
          else
            {
              size = source->diffuseSize;
            }
          unsigned aCoord(source->adjoins[
                        gsl_rng_uniform_int(theRng, size)]);
          Voxel* target(&theLattice[aCoord]);
          if(target->id == theVacantID)
            {
              if(theWalkProbability == 1 ||
                 gsl_rng_uniform(theRng) < theWalkProbability)
                {
                  target->id = theID;
                  source->id = theVacantID;
                  theCoords[i] = aCoord;
                }
            }
        }
    }
  void react(unsigned srcCoord, unsigned tarCoord, unsigned sourceIndex,
             unsigned targetIndex, Species* targetSpecies)
    {
      DiffusionInfluencedReactionProcess* aReaction(
               theDiffusionInfluencedReactions[targetSpecies->getID()]);
      unsigned moleculeA(srcCoord);
      unsigned moleculeB(tarCoord);
      unsigned indexA(sourceIndex);
      unsigned indexB(targetIndex);
      if(aReaction->getA() != this)
        {
          indexA = targetIndex; 
          indexB = sourceIndex;
          moleculeA = tarCoord;
          moleculeB = srcCoord;
        }
      if(aReaction->react(moleculeA, moleculeB, indexA, indexB))
        {
          //Soft remove the source molecule, i.e., keep the id:
          softRemoveMoleculeIndex(sourceIndex);
          //Soft remove the target molecule:
          //Make sure the targetIndex is valid:
          //Target and Source are same species:
          //For some reason if I use theCoords[sourceIndex] instead
          //of getCoord(sourceIndex) the walk method becomes
          //much slower when it is only diffusing without reacting:
          if(targetSpecies == this && getCoord(sourceIndex) == tarCoord)
            {
              softRemoveMoleculeIndex(sourceIndex);
            }
          //If the targetSpecies is a multiscale species with implicit
          //molecule, theTargetIndex is equal to the target molecule size,
          //so we use this info to avoid removing the implicit target molecule:
          else if(targetIndex != targetSpecies->size())
            {
              targetSpecies->softRemoveMoleculeIndex(targetIndex);
            }
          theFinalizeReactions[targetSpecies->getID()] = true;
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
      for(unsigned i(0); i < theMoleculeSize; ++i) 
        {
          unsigned aCoord(getCoord(i));
          if(theStepper->isRemovableEdgeCoord(aCoord, theComp))
            {
              Comp* aSuperComp(
                 theStepper->system2Comp(theComp->system->getSuperSystem())); 
              aSuperComp->vacantSpecies->addCompVoxel(aCoord);
            }
          else 
            { 
              theCoords[newCoordSize] = aCoord;
              ++newCoordSize; 
            }
        }
      theMoleculeSize = newCoordSize;
      //Must resize, otherwise compVoxelSize will be inaccurate:
      theCoords.resize(theMoleculeSize);
      theVariable->setValue(theMoleculeSize);
    }
  void removePeriodicEdgeVoxels()
    {
      int newCoordSize(0);
      for(unsigned i(0); i < theMoleculeSize; ++i) 
        {
          unsigned aCoord(getCoord(i));
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
      theMoleculeSize = newCoordSize;
      //Must resize, otherwise compVoxelSize will be inaccurate:
      theCoords.resize(theMoleculeSize);
      theVariable->setValue(theMoleculeSize);
    }
  void updateSpecies()
    {
      if(isCompVacant && (isDiffusiveVacant || isReactiveVacant))
        {
          theCompCoords = new std::vector<unsigned>;
          for(unsigned i(0); i != theMoleculeSize; ++i)
            { 
              theCompCoords->push_back(theCoords[i]);
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
      else if(isTag)
        {
          updateTagMolecules();
        }
    }
  //If it isReactiveVacant it will only be called by SNRP when it is substrate:
  void updateMoleculeSize()
    {
      if(isDiffusiveVacant || isReactiveVacant)
        {
          updateVacantCoordSize();
        }
    }
  void updateTagMolecules()
    {
      theMoleculeSize = 0;
      for(unsigned i(0); i != theTaggedSpeciesList.size(); ++i)
        {
          Species* aSpecies(theTaggedSpeciesList[i]);
          for(unsigned j(0); j != aSpecies->size(); ++j)
            {
              if(aSpecies->getTagID(j) == theID)
                {
                  unsigned aCoord(aSpecies->getCoord(j));
                  ++theMoleculeSize;
                  if(theMoleculeSize > theCoords.size())
                    {
                      theCoords.push_back(aCoord);
                    }
                  else
                    {
                      theCoords[theMoleculeSize-1] = aCoord;
                    }
                }
            }
        }
    }
  //Even if it is a isCompVacant, this method will be called by
  //VisualizationLogProcess, or SNRP if it is Reactive, or DiffusionProcess
  //if it is Diffusive:
  void updateVacantMolecules()
    {
      theMoleculeSize = 0;
      int aSize(theVacantSpecies->compVoxelSize());
      for(int i(0); i != aSize; ++i)
        { 
          //Voxel* aVoxel(theVacantSpecies->getCompVoxel(i));
          unsigned aCoord(theVacantSpecies->getCompCoord(i));
          //if(aVoxel->id == theID)
          if(theLattice[aCoord].id == theID)
            {
              ++theMoleculeSize;
              if(theMoleculeSize > theCoords.size())
                {
                  theCoords.push_back(aCoord);
                }
              else
                {
                  theCoords[theMoleculeSize-1] = aCoord;
                }
            }
        }
      theVariable->setValue(theMoleculeSize);
    }
  void updateVacantCoordSize()
    {
      theMoleculeSize = 0;
      int aSize(theVacantSpecies->compVoxelSize());
      for(int i(0); i != aSize; ++i)
        { 
          //Voxel* aVoxel(theVacantSpecies->getCompVoxel(i));
          unsigned aCoord(theVacantSpecies->getCompCoord(i));
          //if(aVoxel->id == theID)
          if(theLattice[aCoord].id == theID)
            {
              ++theMoleculeSize;
            }
        }
      if(theMoleculeSize > theCoords.size())
        {
          theCoords.resize(theMoleculeSize);
        }
      theVariable->setValue(theMoleculeSize);
    }
  void setTagID(unsigned anIndex, unsigned anID)
    {
      theTags[anIndex].id = anID;
    }
  unsigned getTagID(unsigned anIndex)
    {
      return theTags[anIndex].id;
    }
  Tag& getTag(unsigned anIndex)
    {
      if(isTagged && anIndex != theMoleculeSize)
        {
          return theTags[anIndex];
        }
      return theNullTag;
    }
  void addMolecule(unsigned aCoord)
    {
      addMolecule(aCoord, theNullTag);
    }
  Species* getMultiscaleVacantSpecies()
    {
      return theMultiscaleVacantSpecies;
    }
  void addMolecule(unsigned aCoord, Tag& aTag)
    {
      if(isMultiscale)
        {
          Species* aSpecies(theStepper->id2species(theLattice[aCoord].id));
          if(aSpecies->getVacantSpecies() != theMultiscaleVacantSpecies)
            {
              doAddMolecule(aCoord, aTag);
              addMultiscaleMolecule(aCoord);
            }
        }
      else if(!isVacant)
        {
          doAddMolecule(aCoord, aTag);
        }
      theLattice[aCoord].id = theID;
    }
  void doAddMolecule(unsigned aCoord, Tag& aTag)
    {
      ++theMoleculeSize; 
      if(theMoleculeSize > theCoords.size())
        {
          theCoords.resize(theMoleculeSize);
          theTags.resize(theMoleculeSize);
        }
      theCoords[theMoleculeSize-1] = aCoord;
      if(isTagged)
        {
          //If it is theNullTag:
          if(aTag.origin == theNullCoord)
            {
              Tag aNewTag = {getCoord(theMoleculeSize-1), theNullID};
              theTags[theMoleculeSize-1] = aNewTag;
            }
          else
            {
              theTags[theMoleculeSize-1] = aTag;
            }
        }
      theVariable->setValue(theMoleculeSize);
    }
  void addMultiscaleMolecule(unsigned aCoord)
    {
      unsigned coordA(aCoord-vacStartCoord);
      for(unsigned i(0); i != theIntersectLipids[coordA].size(); ++i)
        {
          unsigned coordB(theIntersectLipids[coordA][i]+lipStartCoord);
          if(theLattice[coordB].id == theMultiscaleVacantSpecies->getID())
            {
              theLattice[coordB].id = theID;
            }
          else
            {
              multiscaleBind(coordB);
            }
        }
    }
  void removeMultiscaleMolecule(unsigned aCoord)
    {
      unsigned coordA(aCoord-vacStartCoord);
      for(unsigned i(0); i != theIntersectLipids[coordA].size(); ++i)
        {
          unsigned coordB(theIntersectLipids[coordA][i]+lipStartCoord);
          if(theLattice[coordB].id == theID)
            {
              theLattice[coordB].id = theMultiscaleVacantSpecies->getID();
            }
          else
            {
              multiscaleUnbind(coordB);
            }
        }
    }
  bool isIntersectMultiscale(unsigned aCoord)
    {
      unsigned coordA(aCoord-vacStartCoord);
      for(unsigned i(0); i != theIntersectLipids[coordA].size(); ++i)
        {
          unsigned coordB(theIntersectLipids[coordA][i]+lipStartCoord);
          unsigned anID(theLattice[coordB].id);
          if(anID == theID ||
             std::find(theMultiscaleIDs.begin(), theMultiscaleIDs.end(),
                       anID) != theMultiscaleIDs.end())
            {
              return true;
            }
        }
      return false;
    }
  bool isIntersectMultiscale(unsigned srcCoord, unsigned tarCoord)
    {
      bool isIntersect(false);
      unsigned coordA(srcCoord-vacStartCoord);
      std::vector<unsigned> temp;
      temp.resize(theIntersectLipids[coordA].size());
      for(unsigned i(0); i != theIntersectLipids[coordA].size(); ++i)
        {
          unsigned coordB(theIntersectLipids[coordA][i]+lipStartCoord);
          temp[i] = theLattice[coordB].id;
          theLattice[coordB].id = theSpeciesSize;
        }
      isIntersect = isIntersectMultiscale(tarCoord);
      coordA = srcCoord-vacStartCoord;
      for(unsigned i(0); i != theIntersectLipids[coordA].size(); ++i)
        {
          unsigned coordB(theIntersectLipids[coordA][i]+lipStartCoord);
          theLattice[coordB].id = temp[i];
        }
      return isIntersect;
    }
  void multiscaleBind(unsigned aCoord)
    {
      unsigned anID(theLattice[aCoord].id);
      Species* source(theStepper->id2species(anID));
      Species* target(theStepper->id2species(theMultiscaleBindIDs[anID]));
      source->softRemoveMolecule(aCoord);
      target->addMolecule(aCoord);
    }
  void multiscaleUnbind(unsigned aCoord)
    {
      unsigned anID(theLattice[aCoord].id);
      Species* source(theStepper->id2species(anID));
      Species* target(theStepper->id2species(theMultiscaleUnbindIDs[anID]));
      source->softRemoveMolecule(aCoord);
      target->addMolecule(aCoord);
    }
  void addCompVoxel(unsigned aCoord)
    {
      theLattice[aCoord].id = theID;
      //theCompVoxels->push_back(&theLattice[aCoord]);
      theCompCoords->push_back(aCoord);
      ++theMoleculeSize;
      theVariable->setValue(theMoleculeSize);
    }
  String getIDString(Species* aSpecies)
    {
      Variable* aVariable(aSpecies->getVariable());
      if(aVariable)
        {
          return "["+aVariable->getSystemPath().asString()+":"+
            aVariable->getID()+"]["+int2str(aSpecies->getID())+"]";
        }
      else if(aSpecies->getID() == theNullID)
        {
          return "[theNullID]["+int2str(aSpecies->getID())+"]";
        }
      return "[unknown]";
    }
  String getIDString()
    {
      Variable* aVariable(getVariable());
      if(aVariable)
        {
          return "["+aVariable->getSystemPath().asString()+":"+
            aVariable->getID()+"]["+int2str(getID())+"]";
        }
      else if(getID() == theNullID)
        {
          return "[theNullID]["+int2str(getID())+"]";
        }
      return "[unknown]";
    }
  unsigned compVoxelSize()
    {
      //return theCompVoxels->size();
      return theCompCoords->size();
    }
  unsigned getCompCoord(unsigned index)
    {
      return (*theCompCoords)[index];
    }
  /*
  Voxel* getCompVoxel(unsigned index)
    {
      return (*theCompVoxels)[index];
    }
    */
  unsigned getMoleculeIndexFast(unsigned aCoord)
    {
      for(unsigned i(0); i < theMoleculeSize; ++i)
        {
          if(theCoords[i] == aCoord)
            {
              return i;
            }
        }
      return theMoleculeSize;
    }
  unsigned getMoleculeIndex(unsigned aCoord)
    {
      unsigned index(getMoleculeIndexFast(aCoord));
      if(index == theMoleculeSize)
        { 
          if(isDiffusiveVacant || isReactiveVacant)
            {
              updateVacantMolecules();
            }
          index = getMoleculeIndexFast(aCoord);
          if(index == theMoleculeSize)
            { 
              std::cout << "error in getting the index:" << getIDString() <<
               " size:" << theMoleculeSize << std::endl;
              return 0;
            }
        }
      return index;
    }
  //it is soft remove because the id of the molecule is not changed:
  void softRemoveMolecule(unsigned aCoord)
    {
      if(isMultiscale)
        {
          Species* aSpecies(theStepper->id2species(theLattice[aCoord].id));
          if(aSpecies->getVacantSpecies() != theMultiscaleVacantSpecies)
            {
              softRemoveMoleculeIndex(getMoleculeIndex(aCoord));
            }
        }
      else if(!isVacant)
        {
          softRemoveMoleculeIndex(getMoleculeIndex(aCoord));
        }
    }
  void removeMolecule(unsigned aCoord)
    {
      if(isMultiscale)
        {
          Species* aSpecies(theStepper->id2species(theLattice[aCoord].id));
          if(aSpecies->getVacantSpecies() != theMultiscaleVacantSpecies)
            {
              removeMoleculeIndex(getMoleculeIndex(aCoord));
            }
        }
      if(!isVacant)
        {
          removeMoleculeIndex(getMoleculeIndex(aCoord));
        }
    }
  void removeMoleculeIndex(unsigned anIndex)
    {
      if(!isVacant)
        {
          theLattice[theCoords[anIndex]].id = theVacantID;
          softRemoveMoleculeIndex(anIndex);
        }
    }
  void softRemoveMoleculeIndex(unsigned anIndex)
    {
      if(isMultiscale)
        {
          removeMultiscaleMolecule(theCoords[anIndex]);
        }
      if(!isVacant)
        {
          theCoords[anIndex] = theCoords[--theMoleculeSize];
          if(isTagged)
            {
              theTags[anIndex] = theTags[theMoleculeSize];
            }
          theVariable->setValue(theMoleculeSize);
          return;
        }
    }
  //Used to remove all molecules and free memory used to store the molecules
  void clearMolecules()
    {
      theCoords.resize(0);
      theMoleculeSize = 0;
      theVariable->setValue(0);
    }
  //Used by the SpatiocyteStepper when resetting an interation, so must
  //clear the whole compartment using theComp->vacantSpecies->getVacantID():
  void removeMolecules()
    {
      if(!isCompVacant)
        {
          for(unsigned i(0); i < theMoleculeSize; ++i)
            {
              theLattice[theCoords[i]].id = theVacantSpecies->getID();
            }
          theMoleculeSize = 0;
          theVariable->setValue(theMoleculeSize);
        }
    }
  int getPopulateCoordSize()
    {
      return theInitCoordSize-theMoleculeSize;
    }
  int getInitCoordSize()
    {
      return theInitCoordSize;
    }
  void initMoleculeOrigins()
    {
      theMoleculeOrigins.resize(theMoleculeSize);
      for(unsigned i(0); i < theMoleculeSize; ++i)
        {
          Origin& anOrigin(theMoleculeOrigins[i]);
          anOrigin.point = theStepper->coord2point(getCoord(i));
          anOrigin.row = 0;
          anOrigin.layer = 0;
          anOrigin.col = 0;
        }
    }
  void removeBoundaryMolecules()
    {
      for(unsigned i(0); i < theMoleculeSize; ++i)
        {
          if(theStepper->isBoundaryCoord(getCoord(i), theDimension))
            {
              std::cout << "is still there" << std::endl;
            }
        }
      theVariable->setValue(theMoleculeSize);
    }
  void relocateBoundaryMolecules()
    {
      for(unsigned i(0); i < theMoleculeSize; ++i)
        {
          Origin anOrigin(theMoleculeOrigins[i]);
          unsigned periodicCoord(theStepper->getPeriodicCoord(
                                                getCoord(i),
                                                theDimension, &anOrigin));
          if(theLattice[periodicCoord].id == theVacantID)
            {
              theLattice[theCoords[i]].id = theVacantID;
              theCoords[i] = periodicCoord;
              theLattice[theCoords[i]].id = theID;
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
  double getMoleculeRadius() const
    {
      return theMoleculeRadius;
    }
  double getDiffuseRadius() const
    {
      return theDiffuseRadius;
    }
  void setMoleculeRadius(double aRadius)
    {
      theMoleculeRadius = aRadius;
      theDiffuseRadius = aRadius;
    }
  void setDiffuseRadius(double aRadius)
    {
      theDiffuseRadius = aRadius;
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
      if(anInterval < theDiffusionInterval)
        {
          theDiffusionInterval = anInterval;
        }
      for(unsigned i(0); i != theTagSpeciesList.size(); ++i)
        {
          theTagSpeciesList[i]->setDiffusionInterval(theDiffusionInterval);
        }
    }
  unsigned getRandomIndex()
    {
      return gsl_rng_uniform_int(theRng, theMoleculeSize);
    }
  unsigned getRandomMolecule()
    {
      return theCoords[getRandomIndex()];
    }
  void addInterruptedProcess(SpatiocyteNextReactionProcess* aProcess)
    {
      if(std::find(theInterruptedProcesses.begin(),
                   theInterruptedProcesses.end(), aProcess) == 
         theInterruptedProcesses.end())
        {
          theInterruptedProcesses.push_back(aProcess);
        }
    }
  int getBendIndex(double aBendAngle)
    {
      for(unsigned i(0); i != theBendAngles.size(); ++i)
        {
          if(theBendAngles[i] == aBendAngle)
            {
              return i;
            }
        }
      return 0;
    }
  unsigned getRandomAdjoin(unsigned srcCoord, int searchVacant)
    {
      Voxel& source(theLattice[srcCoord]);
      std::vector<unsigned> compCoords;
      if(searchVacant)
        { 
          for(unsigned i(0); i != source.adjoinSize; ++i)
            {
              unsigned aCoord(source.adjoins[i]);
              if(isPopulatable(aCoord))
                {
                  compCoords.push_back(aCoord);
                }
            }
        }
      else
        {
          for(unsigned i(0); i != source.adjoinSize; ++i)
            {
              unsigned aCoord(source.adjoins[i]);
              if(theStepper->id2Comp(theLattice[aCoord].id) == theComp)
                {
                  compCoords.push_back(aCoord);
                }
            }
        }
      return getRandomVacantCoord(compCoords);
    } 
  unsigned getBindingSiteAdjoin(unsigned srcCoord, int bindingSite)
    {
      Voxel& source(theLattice[srcCoord]);
      if(bindingSite < source.adjoinSize)
        { 
          unsigned aCoord(source.adjoins[bindingSite]);
          if(isPopulatable(aCoord))
            {
              return aCoord;
            }
        }
      return theNullCoord;
    } 
  unsigned getRandomAdjoin(unsigned srcCoord, Species* aTargetSpecies,
                           int searchVacant)
    {
      Voxel& source(theLattice[srcCoord]);
      std::vector<unsigned> compCoords;
      if(searchVacant)
        { 
          for(unsigned i(0); i != source.adjoinSize; ++i)
            {
              unsigned aCoord(source.adjoins[i]);
              if(theLattice[aCoord].id == aTargetSpecies->getID())
                {
                  compCoords.push_back(aCoord);
                }
            }
        }
      else
        {
          for(unsigned i(0); i != source.adjoinSize; ++i)
            {
              unsigned aCoord(source.adjoins[i]);
              if(theStepper->id2Comp(theLattice[aCoord].id) == theComp)
                {
                  compCoords.push_back(aCoord);
                }
            }
        }
      return getRandomVacantCoord(compCoords, aTargetSpecies);
    } 
  unsigned getRandomAdjoin(unsigned srcCoord, unsigned tarCoord,
                           int searchVacant)
    {
      Voxel& source(theLattice[srcCoord]);
      std::vector<unsigned> compCoords;
      if(searchVacant)
        { 
          for(unsigned i(0); i != source.adjoinSize; ++i)
            {
              unsigned aCoord(source.adjoins[i]);
              if(aCoord != tarCoord && isPopulatable(aCoord))
                {
                  compCoords.push_back(aCoord);
                }
            }
        }
      else
        {
          for(unsigned i(0); i != source.adjoinSize; ++i)
            {
              unsigned aCoord(source.adjoins[i]);
              if(theStepper->id2Comp(theLattice[aCoord].id) == theComp &&
                 aCoord != tarCoord)
                {
                  compCoords.push_back(aCoord);
                }
            }
        }
      return getRandomVacantCoord(compCoords);
    }
  unsigned getAdjoinMoleculeCnt(unsigned srcCoord, Species* aTargetSpecies)
    {
      Voxel& source(theLattice[srcCoord]);
      unsigned cnt(0);
      for(unsigned i(0); i != source.adjoinSize; ++i)
        {
          if(theLattice[source.adjoins[i]].id == aTargetSpecies->getID())
            {
              ++cnt;
            }
        }
      return cnt;
    }
  unsigned getRandomAdjoin(unsigned srcCoord, unsigned targetA,
                           unsigned targetB, int searchVacant)
    {
      Voxel& source(theLattice[srcCoord]);
      std::vector<unsigned> compCoords;
      if(srcCoord != targetA && srcCoord != targetB)
        {
          if(searchVacant)
            { 
              for(unsigned i(0); i != source.adjoinSize; ++i)
                {
                  unsigned aCoord(source.adjoins[i]);
                  if(isPopulatable(aCoord))
                    {
                      compCoords.push_back(aCoord);
                    }
                }
            }
          else
            {
              for(unsigned i(0); i != source.adjoinSize; ++i)
                {
                  unsigned aCoord(source.adjoins[i]);
                  if(theStepper->id2Comp(theLattice[aCoord].id) == theComp)
                    {
                      compCoords.push_back(aCoord);
                    }
                }
            }
        }
      return getRandomVacantCoord(compCoords);
    }
  unsigned getRandomVacantCoord(std::vector<unsigned>& aCoords)
    {
      if(aCoords.size())
        {
          const int r(gsl_rng_uniform_int(theRng, aCoords.size())); 
          unsigned aCoord(aCoords[r]);
          if(isPopulatable(aCoord))
            {
              return aCoord;
            }
        }
      return theNullCoord;
    }
  unsigned getRandomVacantCoord(std::vector<unsigned>& aCoords,
                                Species* aVacantSpecies)
    {
      if(aCoords.size())
        {
          const int r(gsl_rng_uniform_int(theRng, aCoords.size())); 
          unsigned aCoord(aCoords[r]);
          if(theLattice[aCoord].id == aVacantSpecies->getID())
            {
              return aCoord;
            }
        }
      return theNullCoord;
    }
  unsigned getRandomCompCoord(int searchVacant)
    {
      Species* aVacantSpecies(theComp->vacantSpecies);
      int aSize(aVacantSpecies->compVoxelSize());
      int r(gsl_rng_uniform_int(theRng, aSize));
      if(searchVacant)
        {
          for(int i(r); i != aSize; ++i)
            {
              unsigned aCoord(aVacantSpecies->getCompCoord(i));
              if(isPopulatable(aCoord))
                {
                  return aCoord;
                }
            }
          for(int i(0); i != r; ++i)
            {
              unsigned aCoord(aVacantSpecies->getCompCoord(i));
              if(isPopulatable(aCoord))
                {
                  return aCoord;
                }
            }
        }
      else
        {
          unsigned aCoord(aVacantSpecies->getCompCoord(r));
          if(isPopulatable(aCoord))
            {
              return aCoord;
            }
        }
      return theNullCoord;
    }
  unsigned getRandomAdjoinCompCoord(Comp* aComp, int searchVacant)
    {
      int aSize(theVacantSpecies->size());
      int r(gsl_rng_uniform_int(theRng, aSize)); 
      return getRandomAdjoin(theVacantSpecies->getCoord(r), searchVacant);
    }
  //We need to updateMolecules to set the valid address of voxels
  //since they may have been changed when theLattice is resized by 
  //processes:
  void updateMoleculePointers()
    {
    }
  void saveCoords()
    {
    }
  void setVacStartCoord(unsigned aCoord)
    {
      vacStartCoord = aCoord;
    }
  void setLipStartCoord(unsigned aCoord)
    {
      lipStartCoord = aCoord;
    }
  void setIntersectLipids(Species* aLipid)
    {
      //Traverse through the entire compartment voxels:
      unsigned endA(vacStartCoord+theVacantSpecies->size());
      unsigned endB(lipStartCoord+aLipid->size());
      double dist((aLipid->getMoleculeRadius()+theMoleculeRadius)/
                  (2*theVoxelRadius));
      theIntersectLipids.resize(theVacantSpecies->size());
      for(unsigned i(vacStartCoord); i != endA; ++i)
        {
          Point& pointA(*theLattice[i].point);
          for(unsigned j(lipStartCoord); j != endB; ++j)
            {
              Point& pointB(*theLattice[j].point);
              if(getDistance(pointA, pointB) < dist)
                {
                  //We save j-lipStartCoord and not the absolute coord
                  //since the absolute coord may change after resize 
                  //of lattice:
                  theIntersectLipids[i-vacStartCoord
                    ].push_back(j-lipStartCoord);
                }
            }
        }
    }
  void setMultiscaleBindUnbindIDs(unsigned anID, unsigned aPairID)
    {
      if(std::find(theMultiscaleIDs.begin(), theMultiscaleIDs.end(),
                   anID) == theMultiscaleIDs.end())
        {
          theMultiscaleIDs.push_back(anID);
        }
      theMultiscaleBindIDs[aPairID] = anID;
      theMultiscaleUnbindIDs[anID] = aPairID;
    }
  //Get the fraction of number of nanoscopic molecules (anID) within the
  //multiscale molecule (index):
  double getMultiscaleBoundFraction(unsigned index, unsigned anID)
    {
      double fraction(0);
      if(isMultiscale)
        {
          unsigned i(getCoord(index)-vacStartCoord);
          for(unsigned j(0); j != theIntersectLipids[i].size(); ++j)
            {
              unsigned aCoord(theIntersectLipids[i][j]+lipStartCoord);
              if(theLattice[aCoord].id == anID)
                {
                  fraction += 1;
                }
            }
          fraction /= theIntersectLipids[i].size();
        }
      return fraction;
    }
  unsigned getPopulatableSize()
    {
      if(isMultiscale)
        {
          if(!getIsPopulated())
            {
              std::cout << "The multiscale species:" << 
                getVariable()->getFullID().asString() << " has not yet " <<
                "been populated, but it being populated on." << std::endl;
            }
          thePopulatableCoords.resize(0);
          for(unsigned i(0); i != theMoleculeSize; ++i)
            {
              unsigned j(getCoord(i)-vacStartCoord);
              for(unsigned k(0); k != theIntersectLipids[j].size(); ++k)
                {
                  unsigned aCoord(theIntersectLipids[j][k]+lipStartCoord);
                  if(theLattice[aCoord].id == theID)
                    {
                      thePopulatableCoords.push_back(aCoord);
                    }
                }
            }
          return thePopulatableCoords.size();
        }
      return theMoleculeSize;
    }
  unsigned getRandomPopulatableMolecule()
    {
      unsigned aMolecule;
      if(isMultiscale)
        {
          unsigned index(0);
          do
            {
              index = gsl_rng_uniform_int(theRng, thePopulatableCoords.size());
            }
          while(theLattice[thePopulatableCoords[index]].id != theID);
          aMolecule =  thePopulatableCoords[index];
        }
      else
        {
          aMolecule = getRandomMolecule();
          while(theLattice[aMolecule].id != theID)
            {
              aMolecule = getRandomMolecule();
            }
        }
      return aMolecule;
    }
  unsigned getPopulatableCoord(unsigned index)
    {
      if(isMultiscale)
        {
          return thePopulatableCoords[index];
        }
      return getCoord(index);
    }
  Point coord2point(unsigned aCoord)
    {
      if(theLattice[aCoord].point)
        {
          return *theLattice[aCoord].point;
        }
      return theStepper->coord2point(aCoord);
    }
  void setMultiscaleVacantSpecies(Species* aSpecies)
    {
      theMultiscaleVacantSpecies = aSpecies;
    }
  //Can aVoxel be populated by this species:
  bool isPopulatable(unsigned aCoord)
    {
      if(isMultiscale)
        {
          if(isIntersectMultiscale(aCoord))
            {
              return false;
            }
        }
      else if(theLattice[aCoord].id != theVacantID)
        {
          return false;
        }
      return true;
    }
  //Can aVoxel of this species replaced by aSpecies:
  bool isReplaceable(unsigned aCoord, Species* aSpecies)
    {
      if(getComp() != aSpecies->getComp() &&
         theID != aSpecies->getVacantID())
        {
          return false;
        }
      if(aSpecies->getIsMultiscale())
        {
          if(aSpecies->isIntersectMultiscale(aCoord))
            {
              return false;
            }
        }
      return true;
    }
private:
  bool isCentered;
  bool isCompVacant;
  bool isDiffusing;
  bool isDiffusiveVacant;
  bool isFixedAdjoins;
  bool isGaussianPopulation;
  bool isInContact;
  bool isMultiscale;
  bool isOffLattice;
  bool isPolymer;
  bool isReactiveVacant;
  bool isSubunitInitialized;
  bool isTag;
  bool isTagged;
  bool isVacant;
  const unsigned short theID;
  unsigned lipStartCoord;
  unsigned theAdjoinSize;
  unsigned theCollision;
  unsigned theDimension;
  unsigned theInitCoordSize;
  unsigned theMoleculeSize;
  unsigned theNullCoord;
  unsigned theNullID;
  unsigned theSpeciesSize;
  unsigned vacStartCoord;
  int thePolymerDirectionality;
  unsigned short theVacantID;
  double D;
  double theDiffuseRadius;
  double theDiffusionInterval;
  double theMoleculeRadius;
  double theVoxelRadius;
  double theWalkProbability;
  const gsl_rng* theRng;
  Species* theVacantSpecies;
  Species* theMultiscaleVacantSpecies;
  Comp* theComp;
  MoleculePopulateProcessInterface* thePopulateProcess;
  SpatiocyteStepper* theStepper;
  Variable* theVariable;
  Tag theNullTag;
  std::vector<bool> theFinalizeReactions;
  std::vector<unsigned> collisionCnts;
  std::vector<unsigned> theCoords;
  std::vector<unsigned> theMultiscaleBindIDs;
  std::vector<unsigned> theMultiscaleIDs;
  std::vector<unsigned> theMultiscaleUnbindIDs;
  std::vector<unsigned> thePopulatableCoords;
  std::vector<Tag> theTags;
  std::vector<double> theBendAngles;
  std::vector<double> theReactionProbabilities;
  std::vector<unsigned>* theCompCoords;
  std::vector<Species*> theDiffusionInfluencedReactantPairs;
  std::vector<Species*> theTaggedSpeciesList;
  std::vector<Species*> theTagSpeciesList;
  std::vector<DiffusionInfluencedReactionProcess*> 
    theDiffusionInfluencedReactions;
  std::vector<SpatiocyteNextReactionProcess*> theInterruptedProcesses;
  std::vector<Origin> theMoleculeOrigins;
  std::vector<Voxel>& theLattice;
  std::vector<std::vector<unsigned> > theIntersectLipids;
};


#endif /* __SpatiocyteSpecies_hpp */

