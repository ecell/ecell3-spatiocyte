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
 *    diffuses. In this case, theMolecules list and theMoleculeSize are
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
    isMultiscale(false),
    isOffLattice(false),
    isPeriodic(false),
    isPolymer(false),
    isReactiveVacant(false),
    isRegularLattice(false),
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
    theCompVoxels(&theMolecules),
    theLattice(aLattice) {}
  ~Species() {}
  void initialize(int speciesSize, int anAdjoiningCoordSize,
                  unsigned aNullCoord, unsigned aNullID)
    {
      theAdjoiningCoordSize = anAdjoiningCoordSize;
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
  void setIsPeriodic()
    {
      isPeriodic = true;
    }
  void setIsRegularLattice(unsigned aDiffuseSize)
    {
      isRegularLattice = true;
      theDiffuseSize = aDiffuseSize;
    }
  bool getIsRegularLattice()
    {
      return isRegularLattice;
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
  void populateCompUniformDense(unsigned* voxelIDs, unsigned* aCount)
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
  std::vector<unsigned> getSourceCoords()
    {
      std::vector<unsigned> aCoords;
      /*
      for(unsigned i(0); i != theMoleculeSize; ++i)
        {
          std::vector<unsigned>& 
            aSourceCoords(theMolecules[i]->subunit->sourceCoords);
          for(unsigned j(0); j != aSourceCoords.size(); ++j)
            {
              if(aSourceCoords[j] != theNullCoord)
                {
                  aCoords.push_back(aSourceCoords[j]);
                }
            }
        }
        */
      return aCoords;
    }
  std::vector<unsigned> getTargetCoords()
    {
      std::vector<unsigned> aCoords;
      /*
      for(unsigned i(0); i != theMoleculeSize; ++i)
        {
          std::vector<unsigned>& 
            aTargetCoords(theMolecules[i]->subunit->targetCoords);
          for(unsigned j(0); j != aTargetCoords.size(); ++j)
            {
              if(aTargetCoords[j] != theNullCoord)
                {
                  aCoords.push_back(aTargetCoords[j]);
                }
            }
        }
        */
      return aCoords;
    }
  std::vector<unsigned> getSharedCoords()
    {
      std::vector<unsigned> aCoords;
      /*
      for(unsigned i(0); i != theMoleculeSize; ++i)
        {
          std::vector<unsigned>& 
            aSharedLipids(theMolecules[i]->subunit->sharedLipids);
          for(unsigned j(0); j != aSharedLipids.size(); ++j)
            {
              if(aSharedLipids[j] != theNullCoord)
                {
                  aCoords.push_back(aSharedLipids[j]);
                }
            }
        }
        */
      return aCoords;
    }
  unsigned size() const
    {
      return theMoleculeSize;
    }
  Voxel* getMolecule(int anIndex)
    {
      return theMolecules[anIndex];
    }
  Point getPoint(int anIndex)
    {
      if(isOffLattice)
        {
          if(theMolecules[anIndex]->point)
            {
              return *theMolecules[anIndex]->point;
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
          double aDistance(getDistance(&theMoleculeOrigins[i].point,
                                       &aCurrentPoint));
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
                  std::random_shuffle(theMolecules.begin(), theMolecules.end());
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
      if(isMultiscale)
        {
          for(unsigned i(0); i != theMultiscaleBoundIDs.size(); ++i)
            {
              std::cout << getIDString(theMultiscaleBoundIDs[i]) << std::endl;
            }
          for(unsigned i(0); i != theMultiscaleBindableIDs.size(); ++i)
            {
              std::cout << getIDString(theMultiscaleBindableIDs[i]) << std::endl;
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
          //theInterruptedProcesses are the processes that will always be
          //interrupted at the end of a walk:
          theInterruptedProcesses[i
            ]->substrateValueChanged(theStepper->getCurrentTime());
        }
    }
  void addCollision(Voxel* aVoxel)
    {
      for(unsigned i(0); i < theMoleculeSize; ++i)
        {
          if(aVoxel == theMolecules[i])
            {
              ++collisionCnts[i];
              return;
            }
        }
      std::cout << "error in species add collision" << std::endl;
    }
  void walk()
    {
      theMolecules.resize(theMoleculeSize);
      std::random_shuffle(theMolecules.begin(), theMolecules.end());
      const unsigned beginMoleculeSize(theMoleculeSize);
      unsigned size(theAdjoiningCoordSize);
      for(unsigned i(0); i < beginMoleculeSize && i < theMoleculeSize; ++i)
        {
          Voxel* source(theMolecules[i]);
          if(!isFixedAdjoins)
            {
              size = source->diffuseSize;
            }
          Voxel* target(&theLattice[source->adjoiningCoords[
                        gsl_rng_uniform_int(theRng, size)]]);
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
          else
            {
              if(target->id == theComp->interfaceID)
                {
                  unsigned coord(gsl_rng_uniform_int(theRng, 
                                                     target->adjoiningSize-
                                                     target->diffuseSize));
                  coord = target->adjoiningCoords[coord+target->diffuseSize];
                  target = &theLattice[coord];
                }
              if(theDiffusionInfluencedReactions[target->id])
                {
                  //If it meets the reaction probability:
                  if(gsl_rng_uniform(theRng) < 
                     theReactionProbabilities[target->id])
                    { 
                      Species* targetSpecies(theStepper->id2species(
                                                                  target->id));
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
                          targetIndex = targetSpecies->getIndex(target);
                        }
                      if(theCollision)
                        { 
                          ++collisionCnts[i];
                          targetSpecies->addCollision(target);
                          if(theCollision != 2)
                            {
                              return;
                            }
                        }
                      unsigned aMoleculeSize(theMoleculeSize);
                      react(source, target, i, targetIndex, targetSpecies);
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
  void walkMultiscalePropensity()
    {
      const unsigned beginMoleculeSize(theMoleculeSize);
      for(unsigned i(0); i < beginMoleculeSize && i < theMoleculeSize; ++i)
        {
          Voxel* source(theMolecules[i]);
          const int size(source->diffuseSize);
          Voxel* target(&theLattice[source->adjoiningCoords[
                        gsl_rng_uniform_int(theRng, size)]]);
          if(target->id == theVacantID)
            {
              if(!isIntersectMultiscale(source->coord, target->coord) &&
                 isMultiscaleWalkPropensity(source->coord, target->coord))
                {
                  removeMultiscaleMolecule(source);
                  addMultiscaleMolecule(target);
                  target->id = theID;
                  source->id = theVacantID;
                  theMolecules[i] = target;
                }
            }
        }
    }
  void walkMultiscalePropensityRegular()
    {
      const unsigned beginMoleculeSize(theMoleculeSize);
      for(unsigned i(0); i < beginMoleculeSize && i < theMoleculeSize; ++i)
        {
          Voxel* source(theMolecules[i]);
          const unsigned tarIndex(gsl_rng_uniform_int(theRng, theDiffuseSize));
          Voxel* target(&theLattice[source->adjoiningCoords[tarIndex]]);
          if(target->id == theVacantID)
            {
              if(!isIntersectMultiscaleRegular(source->coord, tarIndex) &&
                 isMultiscaleWalkPropensityRegular(source->coord, tarIndex))
                {
                  moveMultiscaleMoleculeRegular(source->coord, tarIndex);
                  target->id = theID;
                  source->id = theVacantID;
                  theMolecules[i] = target;
                }
            }
        }
    }
  void walkMultiscale()
    {
      const unsigned beginMoleculeSize(theMoleculeSize);
      for(unsigned i(0); i < beginMoleculeSize && i < theMoleculeSize; ++i)
        {
          Voxel* source(theMolecules[i]);
          const int size(source->diffuseSize);
          Voxel* target(&theLattice[source->adjoiningCoords[
                        gsl_rng_uniform_int(theRng, size)]]);
          if(target->id == theVacantID)
            {
              if(!isIntersectMultiscale(source->coord, target->coord))
                {
                  removeMultiscaleMolecule(source);
                  addMultiscaleMolecule(target);
                  target->id = theID;
                  source->id = theVacantID;
                  theMolecules[i] = target;
                }
            }
        }
    }
  void walkMultiscaleRegular()
    {
      const unsigned beginMoleculeSize(theMoleculeSize);
      for(unsigned i(0); i < beginMoleculeSize && i < theMoleculeSize; ++i)
        {
          Voxel* source(theMolecules[i]);
          const unsigned tarIndex(gsl_rng_uniform_int(theRng, theDiffuseSize));
          Voxel* target(&theLattice[source->adjoiningCoords[tarIndex]]);
          if(target->id == theVacantID)
            {
              if(!isIntersectMultiscaleRegular(source->coord, tarIndex))
                {
                  moveMultiscaleMoleculeRegular(source->coord, tarIndex);
                  target->id = theID;
                  source->id = theVacantID;
                  theMolecules[i] = target;
                }
            }
        }
    }
  void walkVacant()
    {
      updateVacantMolecules();
      for(unsigned i(0); i < theMoleculeSize; ++i)
        {
          Voxel* source(theMolecules[i]);
          int size;
          if(isFixedAdjoins)
            {
              size = theAdjoiningCoordSize;
            }
          else
            {
              size = source->diffuseSize;
            }
          Voxel* target(&theLattice[source->adjoiningCoords[
                        gsl_rng_uniform_int(theRng, size)]]);
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
        }
    }
  void react(Voxel* source, Voxel* target, unsigned sourceIndex,
             unsigned targetIndex, Species* targetSpecies)
    {
      DiffusionInfluencedReactionProcess* aReaction(
               theDiffusionInfluencedReactions[targetSpecies->getID()]);
      Voxel* moleculeA(source);
      Voxel* moleculeB(target);
      unsigned indexA(sourceIndex);
      unsigned indexB(targetIndex);
      if(aReaction->getA() != this)
        {
          indexA = targetIndex; 
          indexB = sourceIndex;
          moleculeA = target;
          moleculeB = source;
        }
      if(aReaction->react(moleculeA, moleculeB, indexA, indexB))
        {
          //Soft remove the source molecule, i.e., keep the id:
          softRemoveMolecule(sourceIndex);
          //Soft remove the target molecule:
          //Make sure the targetIndex is valid:
          //Target and Source are same species:
          //For some reason if I use theMolecules[sourceIndex] instead
          //of getMolecule(sourceIndex) the walk method becomes
          //much slower when it is only diffusing without reacting:
          if(targetSpecies == this && getMolecule(sourceIndex) == target)
            {
              softRemoveMolecule(sourceIndex);
            }
          //If the targetSpecies is a multiscale species with implicit
          //molecule, theTargetIndex is equal to the target molecule size,
          //so we use this info to avoid removing the implicit target molecule:
          else if(targetIndex != targetSpecies->size())
            {
              targetSpecies->softRemoveMolecule(targetIndex);
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
  unsigned getCoord(unsigned anIndex)
    {
      return theMolecules[anIndex]->coord;
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
              theMolecules[newCoordSize] = &theLattice[aCoord];
              ++newCoordSize; 
            }
        }
      theMoleculeSize = newCoordSize;
      //Must resize, otherwise compVoxelSize will be inaccurate:
      theMolecules.resize(theMoleculeSize);
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
              theMolecules[newCoordSize] = &theLattice[aCoord];
              ++newCoordSize; 
            }
        }
      theMoleculeSize = newCoordSize;
      //Must resize, otherwise compVoxelSize will be inaccurate:
      theMolecules.resize(theMoleculeSize);
      theVariable->setValue(theMoleculeSize);
    }
  void updateSpecies()
    {
      if(isCompVacant && (isDiffusiveVacant || isReactiveVacant))
        {
          theCompVoxels = new std::vector<Voxel*>;
          for(unsigned i(0); i != theMoleculeSize; ++i)
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
                  Voxel* aVoxel(aSpecies->getMolecule(j));
                  ++theMoleculeSize;
                  if(theMoleculeSize > theMolecules.size())
                    {
                      theMolecules.push_back(aVoxel);
                    }
                  else
                    {
                      theMolecules[theMoleculeSize-1] = aVoxel;
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
          Voxel* aVoxel(theVacantSpecies->getCompVoxel(i));
          if(aVoxel->id == theID)
            {
              ++theMoleculeSize;
              if(theMoleculeSize > theMolecules.size())
                {
                  theMolecules.push_back(aVoxel);
                }
              else
                {
                  theMolecules[theMoleculeSize-1] = aVoxel;
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
          Voxel* aVoxel(theVacantSpecies->getCompVoxel(i));
          if(aVoxel->id == theID)
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
  void addMolecule(Voxel* aVoxel)
    {
      addMolecule(aVoxel, theNullTag);
    }
  Species* getMultiscaleVacantSpecies()
    {
      return theMultiscaleVacantSpecies;
    }
  void addMolecule(Voxel* aVoxel, Tag& aTag)
    {
      if(isMultiscale)
        {
          Species* aSpecies(theStepper->id2species(aVoxel->id));
          if(aSpecies->getVacantSpecies() != theMultiscaleVacantSpecies)
            {
              doAddMolecule(aVoxel, aTag);
              addMultiscaleMolecule(aVoxel);
            }
        }
      else if(!isVacant)
        {
          doAddMolecule(aVoxel, aTag);
        }
      aVoxel->id = theID;
    }
  void doAddMolecule(Voxel* aVoxel, Tag& aTag)
    {
      ++theMoleculeSize; 
      if(theMoleculeSize > theMolecules.size())
        {
          theMolecules.resize(theMoleculeSize);
          theTags.resize(theMoleculeSize);
        }
      theMolecules[theMoleculeSize-1] = aVoxel;
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
  bool isMultiscaleWalkPropensityRegular(const unsigned srcCoord,
                                         const unsigned tarIndex)
    {
      unsigned tarCnt(0);
      unsigned srcCnt(0);
      const unsigned coordA(srcCoord-vacStartCoord);
      const int rowA(coordA/lipCols);
      const std::vector<int>& anOffsetsA(theTarOffsets[rowA%2][tarIndex]);
      //count tar
      for(unsigned i(0); i != anOffsetsA.size(); ++i)
        {
          const int offsetRow((anOffsetsA[i]+theRegLatticeCoord)/lipCols-
                              theRegLatticeCoord/lipCols);
          int coordB(coordA+anOffsetsA[i]);
          if(isInLattice(coordB, offsetRow+rowA))
            {
              const unsigned anID(theLattice[coordB+lipStartCoord].id);
              if(std::find(theMultiscaleBindableIDs.begin(), 
                           theMultiscaleBindableIDs.end(), anID) !=
                 theMultiscaleBindableIDs.end())
                {
                  ++tarCnt;
                }
            }
        }
      //count src
      const std::vector<int>& anOffsetsB(theSrcOffsets[rowA%2][tarIndex]);
      for(unsigned i(0); i != anOffsetsB.size(); ++i)
        {
          const int offsetRow((anOffsetsB[i]+theRegLatticeCoord)/lipCols-
                              theRegLatticeCoord/lipCols);
          int coordB(coordA+anOffsetsB[i]);
          if(isInLattice(coordB, offsetRow+rowA))
            {
              const unsigned anID(theLattice[coordB+lipStartCoord].id);
              if(std::find(theMultiscaleBoundIDs.begin(), 
                           theMultiscaleBoundIDs.end(), anID) !=
                 theMultiscaleBoundIDs.end())
                {
                  ++srcCnt;
                }
            }
        }
      if(tarCnt > srcCnt)
        {
          return true;
        }
      return false;
    }
  bool isMultiscaleWalkPropensity(const unsigned srcCoord,
                                  const unsigned tarCoord)
    {
      unsigned srcCnt(0);
      unsigned tarCnt(0);
      const unsigned coordA(srcCoord-vacStartCoord);
      const unsigned coordB(tarCoord-vacStartCoord);
      for(unsigned i(0); i != theIntersectLipids[coordA].size(); ++i)
        {
          unsigned coord(theIntersectLipids[coordA][i]+lipStartCoord);
          unsigned anID(theLattice[coord].id);
          if(std::find(theMultiscaleBoundIDs.begin(), 
                       theMultiscaleBoundIDs.end(), anID) !=
             theMultiscaleBoundIDs.end())
            {
              ++srcCnt;
            }
        }
      for(unsigned i(0); i != theIntersectLipids[coordB].size(); ++i)
        {
          unsigned coord(theIntersectLipids[coordB][i]+lipStartCoord);
          unsigned anID(theLattice[coord].id);
          if(std::find(theMultiscaleBoundIDs.begin(), 
                       theMultiscaleBoundIDs.end(), anID) !=
             theMultiscaleBoundIDs.end() ||
             std::find(theMultiscaleBindableIDs.begin(), 
                       theMultiscaleBindableIDs.end(), anID) !=
             theMultiscaleBindableIDs.end())
            {
              ++tarCnt;
            }
        }
      if(tarCnt > srcCnt)
        {
          return true;
        }
      return false;
    }
  void addMultiscaleMolecule(Voxel* aVoxel)
    {
      const unsigned coordA(aVoxel->coord-vacStartCoord);
      if(isRegularLattice)
        {
          const int rowA(coordA/lipCols);
          const std::vector<int>& anOffsets(theOffsets[rowA%2]);
          for(unsigned i(0); i != anOffsets.size(); ++i)
            {
              const int offsetRow((anOffsets[i]+theRegLatticeCoord)/lipCols-
                                  theRegLatticeCoord/lipCols);
              int coordB(coordA+anOffsets[i]);
              if(isInLattice(coordB, offsetRow+rowA))
                {
                  const unsigned coord(coordB+lipStartCoord);
                  if(theLattice[coord].id == 
                     theMultiscaleVacantSpecies->getID())
                    {
                      theLattice[coord].id = theID;
                    }
                  else
                    {
                      multiscaleBind(&theLattice[coord]);
                    }
                }
            }
        }
      else
        {
          for(unsigned i(0); i != theIntersectLipids[coordA].size(); ++i)
            {
              unsigned coordB(theIntersectLipids[coordA][i]+lipStartCoord);
              if(theLattice[coordB].id == theMultiscaleVacantSpecies->getID())
                {
                  theLattice[coordB].id = theID;
                }
              else
                {
                  multiscaleBind(&theLattice[coordB]);
                }
            }
        }
    }
  void removeMultiscaleMolecule(Voxel* aVoxel)
    {
      const unsigned coordA(aVoxel->coord-vacStartCoord);
      if(isRegularLattice)
        {
          const int rowA(coordA/lipCols);
          const std::vector<int>& anOffsets(theOffsets[rowA%2]);
          for(unsigned i(0); i != anOffsets.size(); ++i)
            {
              const int offsetRow((anOffsets[i]+theRegLatticeCoord)/lipCols-
                                  theRegLatticeCoord/lipCols);
              int coordB(coordA+anOffsets[i]);
              if(isInLattice(coordB, offsetRow+rowA))
                {
                  const unsigned coord(coordB+lipStartCoord);
                  if(theLattice[coord].id == theID)
                    {
                      theLattice[coord].id = 
                        theMultiscaleVacantSpecies->getID();
                    }
                  else
                    {
                      multiscaleUnbind(&theLattice[coord]);
                    }
                }
            }
        }
      else
        {
          for(unsigned i(0); i != theIntersectLipids[coordA].size(); ++i)
            {
              unsigned coordB(theIntersectLipids[coordA][i]+lipStartCoord);
              if(theLattice[coordB].id == theID)
                {
                  theLattice[coordB].id = theMultiscaleVacantSpecies->getID();
                }
              else
                {
                  multiscaleUnbind(&theLattice[coordB]);
                }
            }
        }
    }
  bool isIntersectMultiscale(const unsigned srcCoord)
    {
      const unsigned coordA(srcCoord-vacStartCoord);
      if(isRegularLattice)
        {
          const int rowA(coordA/lipCols);
          const std::vector<int>& anOffsets(theOffsets[rowA%2]);
          for(unsigned i(0); i != anOffsets.size(); ++i)
            {
              const int offsetRow((anOffsets[i]+theRegLatticeCoord)/lipCols-
                                  theRegLatticeCoord/lipCols);
              int coordB(coordA+anOffsets[i]);
              if(isInLattice(coordB, offsetRow+rowA))
                {
                  const unsigned anID(theLattice[coordB+lipStartCoord].id);
                  if(anID == theID ||
                     std::find(theMultiscaleBoundIDs.begin(), 
                               theMultiscaleBoundIDs.end(),
                               anID) != theMultiscaleBoundIDs.end())
                    {
                      return true;
                    }
                }
            }
        }
      else
        {
          for(unsigned i(0); i != theIntersectLipids[coordA].size(); ++i)
            {
              unsigned coordB(theIntersectLipids[coordA][i]+lipStartCoord);
              unsigned anID(theLattice[coordB].id);
              if(anID == theID ||
                 std::find(theMultiscaleBoundIDs.begin(), 
                           theMultiscaleBoundIDs.end(),
                           anID) != theMultiscaleBoundIDs.end())
                {
                  return true;
                }
            }
        }
      return false;
    }
  bool isInLattice(int& coord, const int offset)
    {
      if(isPeriodic)
        {
          if(coord < 0)
            {
              coord += lipCols*(offset-(coord+1)/lipCols+1);
            }
          else
            {
              coord += lipCols*(offset-(coord/lipCols));
            }
          if(coord < 0)
            {
              coord += lipRows*lipCols;
            }
          else if(coord >= lipRows*lipCols)
            {
              coord -= lipRows*lipCols;
            }
        }
      else if(coord/lipCols != offset || coord < 0 || coord >= lipRows*lipCols)
        {
          return false;
        }
      return true;
    }
  void moveMultiscaleMoleculeRegular(const unsigned srcCoord, 
                                     const unsigned tarIndex)
    {
      const unsigned coordA(srcCoord-vacStartCoord);
      const int rowA(coordA/lipCols);
      const std::vector<int>& anOffsetsA(theTarOffsets[rowA%2][tarIndex]);
      //Add tar
      for(unsigned i(0); i != anOffsetsA.size(); ++i)
        {
          const int offsetRow((anOffsetsA[i]+theRegLatticeCoord)/lipCols-
                              theRegLatticeCoord/lipCols);
          int coordB(coordA+anOffsetsA[i]);
          if(isInLattice(coordB, offsetRow+rowA))
            {
              const unsigned coord(coordB+lipStartCoord);
              if(theLattice[coord].id == 
                 theMultiscaleVacantSpecies->getID())
                {
                  theLattice[coord].id = theID;
                }
              else
                {
                  multiscaleBind(&theLattice[coord]);
                }
            }
        }
      //Remove src
      const std::vector<int>& anOffsetsB(theSrcOffsets[rowA%2][tarIndex]);
      for(unsigned i(0); i != anOffsetsB.size(); ++i)
        {
          const int offsetRow((anOffsetsB[i]+theRegLatticeCoord)/lipCols-
                              theRegLatticeCoord/lipCols);
          int coordB(coordA+anOffsetsB[i]);
          if(isInLattice(coordB, offsetRow+rowA))
            {
              const unsigned coord(coordB+lipStartCoord);
              if(theLattice[coord].id == theID)
                {
                  theLattice[coord].id = theMultiscaleVacantSpecies->getID();
                }
              else
                {
                  multiscaleUnbind(&theLattice[coord]);
                }
            }
        }
    }
  bool isIntersectMultiscaleRegular(const unsigned srcCoord, 
                                    const unsigned tarIndex)
    {
      const unsigned coordA(srcCoord-vacStartCoord);
      const int rowA(coordA/lipCols);
      const std::vector<int>& anOffsets(theTarOffsets[rowA%2][tarIndex]);
      for(unsigned i(0); i != anOffsets.size(); ++i)
        {
          const int offsetRow((anOffsets[i]+theRegLatticeCoord)/lipCols-
                              theRegLatticeCoord/lipCols);
          int coordB(coordA+anOffsets[i]);
          if(isInLattice(coordB, offsetRow+rowA))
            {
              const unsigned anID(theLattice[coordB+lipStartCoord].id);
              if(anID == theID ||
                 std::find(theMultiscaleBoundIDs.begin(), 
                           theMultiscaleBoundIDs.end(),
                           anID) != theMultiscaleBoundIDs.end())
                {
                  return true;
                }
            }
        }
      return false;
    }
  bool isIntersectMultiscale(const unsigned srcCoord, const unsigned tarCoord)
    {
      bool isIntersect(false);
      const unsigned coordA(srcCoord-vacStartCoord);
      std::vector<unsigned> temp;
      temp.resize(theIntersectLipids[coordA].size());
      for(unsigned i(0); i != theIntersectLipids[coordA].size(); ++i)
        {
          unsigned coordB(theIntersectLipids[coordA][i]+lipStartCoord);
          temp[i] = theLattice[coordB].id;
          theLattice[coordB].id = theSpeciesSize;
        }
      isIntersect = isIntersectMultiscale(tarCoord);
      for(unsigned i(0); i != theIntersectLipids[coordA].size(); ++i)
        {
          unsigned coordB(theIntersectLipids[coordA][i]+lipStartCoord);
          theLattice[coordB].id = temp[i];
        }
      return isIntersect;
    }
  void multiscaleBind(Voxel* aVoxel)
    {
      unsigned anID(aVoxel->id);
      Species* source(theStepper->id2species(anID));
      Species* target(theStepper->id2species(theMultiscaleBindIDs[anID]));
      source->softRemoveMolecule(aVoxel);
      target->addMolecule(aVoxel);
    }
  void multiscaleUnbind(Voxel* aVoxel)
    {
      unsigned anID(aVoxel->id);
      Species* source(theStepper->id2species(anID));
      Species* target(theStepper->id2species(theMultiscaleUnbindIDs[anID]));
      source->softRemoveMolecule(aVoxel);
      target->addMolecule(aVoxel);
    }
  void addCompVoxel(unsigned aCoord)
    {
      theLattice[aCoord].id = theID;
      theCompVoxels->push_back(&theLattice[aCoord]);
      ++theMoleculeSize;
      theVariable->setValue(theMoleculeSize);
    }
  String getIDString(unsigned anID)
    {
      Species* aSpecies(theStepper->id2species(anID));
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
      return theCompVoxels->size();
    }
  Voxel* getCompVoxel(unsigned index)
    {
      return (*theCompVoxels)[index];
    }
  unsigned getIndexFast(Voxel* aVoxel)
    {
      for(unsigned i(0); i < theMoleculeSize; ++i)
        {
          if(theMolecules[i] == aVoxel)
            {
              return i;
            }
        }
      return theMoleculeSize;
    }
  unsigned getIndex(Voxel* aVoxel)
    {
      //This is required by SNRP reactABC with B a lipid or vacant molecule:
      if(getIsCompVacant())
        {
          return theMoleculeSize;
        }
      unsigned index(getIndexFast(aVoxel));
      if(index == theMoleculeSize)
        { 
          if(isDiffusiveVacant || isReactiveVacant)
            {
              updateVacantMolecules();
            }
          index = getIndexFast(aVoxel);
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
  void softRemoveMolecule(Voxel* aVoxel)
    {
      if(isMultiscale)
        {
          Species* aSpecies(theStepper->id2species(aVoxel->id));
          if(aSpecies->getVacantSpecies() != theMultiscaleVacantSpecies)
            {
              softRemoveMolecule(getIndex(aVoxel));
            }
        }
      else if(!isVacant)
        {
          softRemoveMolecule(getIndex(aVoxel));
        }
    }
  void removeMolecule(Voxel* aVoxel)
    {
      if(isMultiscale)
        {
          Species* aSpecies(theStepper->id2species(aVoxel->id));
          if(aSpecies->getVacantSpecies() != theMultiscaleVacantSpecies)
            {
              removeMolecule(getIndex(aVoxel));
            }
        }
      if(!isVacant)
        {
          removeMolecule(getIndex(aVoxel));
        }
    }
  void removeMolecule(unsigned anIndex)
    {
      if(!isVacant)
        {
          theMolecules[anIndex]->id = theVacantID;
          softRemoveMolecule(anIndex);
        }
    }
  void softRemoveMolecule(unsigned anIndex)
    {
      if(isMultiscale)
        {
          removeMultiscaleMolecule(theMolecules[anIndex]);
        }
      if(!isVacant)
        {
          theMolecules[anIndex] = theMolecules[--theMoleculeSize];
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
      theMolecules.resize(0);
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
              theMolecules[i]->id = theVacantSpecies->getID();
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
              theMolecules[i]->id = theVacantID;
              theMolecules[i] = &theLattice[periodicCoord];
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
  Voxel* getRandomMolecule()
    {
      return theMolecules[getRandomIndex()];
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
  Voxel* getRandomAdjoiningVoxel(Voxel* source, int searchVacant)
    {
      std::vector<unsigned> compCoords;
      if(searchVacant)
        { 
          for(unsigned i(0); i != source->adjoiningSize; ++i)
            {
              unsigned aCoord(source->adjoiningCoords[i]);
              if(isPopulatable(&theLattice[aCoord]))
                {
                  compCoords.push_back(aCoord);
                }
            }
        }
      else
        {
          for(unsigned i(0); i != source->adjoiningSize; ++i)
            {
              unsigned aCoord(source->adjoiningCoords[i]);
              if(theStepper->id2Comp(theLattice[aCoord].id) == theComp)
                {
                  compCoords.push_back(aCoord);
                }
            }
        }
      return getRandomVacantVoxel(compCoords);
    } 
  Voxel* getBindingSiteAdjoiningVoxel(Voxel* source, int bindingSite)
    {
      if(bindingSite < source->adjoiningSize)
        { 
          unsigned aCoord(source->adjoiningCoords[bindingSite]);
          if(isPopulatable(&theLattice[aCoord]))
            {
              return &theLattice[aCoord];
            }
        }
      return NULL;
    } 
  Voxel* getRandomAdjoiningVoxel(Voxel* source,
                                 Species* aTargetSpecies,
                                 int searchVacant)
    {
      std::vector<unsigned> compCoords;
      if(searchVacant)
        { 
          for(unsigned i(0); i != source->adjoiningSize; ++i)
            {
              unsigned aCoord(source->adjoiningCoords[i]);
              if(theLattice[aCoord].id == aTargetSpecies->getID())
                {
                  compCoords.push_back(aCoord);
                }
            }
        }
      else
        {
          for(unsigned i(0); i != source->adjoiningSize; ++i)
            {
              unsigned aCoord(source->adjoiningCoords[i]);
              if(theStepper->id2Comp(theLattice[aCoord].id) == theComp)
                {
                  compCoords.push_back(aCoord);
                }
            }
        }
      return getRandomVacantVoxel(compCoords, aTargetSpecies);
    } 
  Voxel* getRandomAdjoiningVoxel(Voxel* source, Voxel* target,
                                 int searchVacant)
    {
      std::vector<unsigned> compCoords;
      if(searchVacant)
        { 
          for(unsigned i(0); i != source->adjoiningSize; ++i)
            {
              unsigned aCoord(source->adjoiningCoords[i]);
              if(&theLattice[aCoord] != target && 
                 isPopulatable(&theLattice[aCoord]))
                {
                  compCoords.push_back(aCoord);
                }
            }
        }
      else
        {
          for(unsigned i(0); i != source->adjoiningSize; ++i)
            {
              unsigned aCoord(source->adjoiningCoords[i]);
              if(theStepper->id2Comp(theLattice[aCoord].id) == theComp &&
                 &theLattice[aCoord] != target)
                {
                  compCoords.push_back(aCoord);
                }
            }
        }
      return getRandomVacantVoxel(compCoords);
    }
  bool isAdjoinedSpecies(Voxel* source, Species* aTargetSpecies)
    {
      for(unsigned i(0); i != source->adjoiningSize; ++i)
        {
          if(theLattice[source->adjoiningCoords[i]].id == 
             aTargetSpecies->getID())
            {
              return true;
            }
        }
      return false;
    }
  unsigned getAdjoiningMoleculeCnt(Voxel* source, Species* aTargetSpecies)
    {
      unsigned cnt(0);
      for(unsigned i(0); i != source->adjoiningSize; ++i)
        {
          if(theLattice[source->adjoiningCoords[i]].id == 
             aTargetSpecies->getID())
            {
              ++cnt;
            }
        }
      return cnt;
    }
  Voxel* getRandomAdjoiningVoxel(Voxel* source, Voxel* targetA, Voxel* targetB,
                                 int searchVacant)
    {
      std::vector<unsigned> compCoords;
      if(searchVacant)
        { 
          for(unsigned i(0); i != source->adjoiningSize; ++i)
            {
              unsigned aCoord(source->adjoiningCoords[i]);
              if(source != targetA && source != targetB && 
                 isPopulatable(&theLattice[aCoord]))
                {
                  compCoords.push_back(aCoord);
                }
            }
        }
      else
        {
          for(unsigned i(0); i != source->adjoiningSize; ++i)
            {
              unsigned aCoord(source->adjoiningCoords[i]);
              if(theStepper->id2Comp(theLattice[aCoord].id) == theComp &&
                 source != targetA && source != targetB)
                {
                  compCoords.push_back(aCoord);
                }
            }
        }
      return getRandomVacantVoxel(compCoords);
    }
  Voxel* getRandomVacantVoxel(std::vector<unsigned>& aCoords)
    {
      if(aCoords.size())
        {
          const int r(gsl_rng_uniform_int(theRng, aCoords.size())); 
          unsigned aCoord(aCoords[r]);
          if(isPopulatable(&theLattice[aCoord]))
            {
              return &theLattice[aCoord];
            }
        }
      return NULL;
    }
  Voxel* getRandomVacantVoxel(std::vector<unsigned>& aCoords,
                              Species* aVacantSpecies)
    {
      if(aCoords.size())
        {
          const int r(gsl_rng_uniform_int(theRng, aCoords.size())); 
          unsigned aCoord(aCoords[r]);
          if(theLattice[aCoord].id == aVacantSpecies->getID())
            {
              return &theLattice[aCoord];
            }
        }
      return NULL;
    }
  Voxel* getRandomCompVoxel(int searchVacant)
    {
      Species* aVacantSpecies(theComp->vacantSpecies);
      int aSize(aVacantSpecies->compVoxelSize());
      int r(gsl_rng_uniform_int(theRng, aSize));
      if(searchVacant)
        {
          for(int i(r); i != aSize; ++i)
            {
              Voxel* aVoxel(aVacantSpecies->getCompVoxel(i));
              if(isPopulatable(aVoxel))
                {
                  return aVoxel;
                }
            }
          for(int i(0); i != r; ++i)
            {
              Voxel* aVoxel(aVacantSpecies->getCompVoxel(i));
              if(isPopulatable(aVoxel))
                {
                  return aVoxel;
                }
            }
        }
      else
        {
          Voxel* aVoxel(aVacantSpecies->getCompVoxel(r));
          if(isPopulatable(aVoxel))
            {
              return aVoxel;
            }
        }
      return NULL;
    }
  Voxel* getRandomAdjoiningCompVoxel(Comp* aComp, int searchVacant)
    {
      int aSize(theVacantSpecies->size());
      int r(gsl_rng_uniform_int(theRng, aSize)); 
      Voxel* aVoxel(theVacantSpecies->getMolecule(r));
      return getRandomAdjoiningVoxel(aVoxel, searchVacant);
    }
  //We need to updateMolecules to set the valid address of voxels
  //since they may have been changed when theLattice is resized by 
  //processes:
  void updateMoleculePointers()
    {
      for(unsigned i(0); i != theCoords.size(); ++i)
        {
          theMolecules[i] = &theLattice[theCoords[i]];
        }
      theCoords.resize(0);
    }
  void saveCoords()
    {
      theCoords.resize(theMoleculeSize);
      for(unsigned i(0); i != theMoleculeSize; ++i)
        {
          theCoords[i] = getCoord(i);
        }
    }
  void setVacStartCoord(unsigned aCoord, unsigned aVacantRows,
                        unsigned aVacantCols)
    {
      vacStartCoord = aCoord;
      vacRows = aVacantRows;
      vacCols = aVacantCols;
    }
  void setLipStartCoord(unsigned aCoord, unsigned aLipidRows,
                        unsigned aLipidCols)
    {
      lipStartCoord = aCoord;
      lipRows = aLipidRows;
      lipCols = aLipidCols;
    }
  void setIntersectLipids(Species* aLipid, Point& aLipidStart, double aGridSize,
                          unsigned aGridCols, unsigned aGridRows, 
                          std::vector<std::vector<unsigned> >& aGrid,
                          unsigned aVacantRows, unsigned aVacantCols)
    {
      double nDist((aLipid->getMoleculeRadius()+theMoleculeRadius)/
                   (2*theVoxelRadius));
      //Traverse through the entire compartment voxels:
      unsigned endA(vacStartCoord+theVacantSpecies->size());
      theIntersectLipids.resize(theVacantSpecies->size());
      for(unsigned i(vacStartCoord); i != endA; ++i)
        {
          getIntersectLipids(i, nDist, aLipidStart, aGridSize,
                             aGridCols, aGridRows, aGrid,
                             theIntersectLipids[i-vacStartCoord]);
        }
    }
  void getIntersectLipids(unsigned coordA, double nDist, Point& aLipidStart,
                          double aGridSize, unsigned aGridCols, 
                          unsigned aGridRows, 
                          std::vector<std::vector<unsigned> >& aGrid,
                          std::vector<unsigned>& anIntersectLipids)
    {
      Point& pointA(*theLattice[coordA].point);
      unsigned row((unsigned)((pointA.y-aLipidStart.y)/aGridSize));
      unsigned col((unsigned)((pointA.z-aLipidStart.z)/aGridSize));
      unsigned rowStart(std::max(unsigned(1), row)-1);
      unsigned rowEnd(std::min(unsigned(aGridRows), row+2));
      for(unsigned j(rowStart); j != rowEnd; ++j)
        {
          unsigned colStart(std::max(unsigned(1), col)-1);
          unsigned colEnd(std::min(unsigned(aGridCols), col+2));
          for(unsigned k(colStart); k != colEnd; ++k)
            {
              std::vector<unsigned>& coords(aGrid[k+aGridCols*j]);
              for(unsigned l(0); l != coords.size(); ++l)
                {
                  unsigned m(coords[l]);
                  Point& pointB(*theLattice[m].point);
                  if(getDistance(&pointA, &pointB) < nDist)
                    {
                      anIntersectLipids.push_back(m-lipStartCoord);
                    }
                }
            }
        }
    }
  void setIntersectOffsets(Species* aLipid, const Point& aLipidStart, 
                                 const unsigned aVacantRows,
                                 const unsigned aVacantCols,
                                 const double nLipidRadius)
    {
      double nDist((aLipid->getMoleculeRadius()+theMoleculeRadius)/
                   (2*theVoxelRadius));
      theRegLatticeCoord = lipRows/2*lipCols+lipCols/2;
      theOffsets.resize(2);
      unsigned rowA(aVacantRows/2);
      unsigned rowB(rowA+1);
      if(rowB >= aVacantRows && rowA > 0)
        {
          rowB = rowA-1;
        }
      else if(!rowA)
        {
          rowB = 0;
        }
      unsigned aCol(aVacantCols/2);
      unsigned coordA(rowA*aVacantCols+aCol);
      unsigned coordB(rowB*aVacantCols+aCol);
      setCoordOffsets(coordA, coordA, nDist, aLipidStart, nLipidRadius,
                      theOffsets[rowA%2]);
      setCoordOffsets(coordB, coordB, nDist, aLipidStart, nLipidRadius,
                      theOffsets[rowB%2]);
      theTarOffsets.resize(2);
      theSrcOffsets.resize(2);
      setWalkOffsets(rowA, coordA, nDist, aLipidStart, nLipidRadius,
                     theOffsets[rowA%2]);
      setWalkOffsets(rowB, coordB, nDist, aLipidStart, nLipidRadius,
                     theOffsets[rowB%2]);
      /*
      std::cout << "row0:" << theDiffuseSize << std::endl;
      for(unsigned i(0); i != theDiffuseSize; ++i)
        {
          std::cout << "tar:" << i << std::endl; 
          for(unsigned j(0); j != theTarOffsets[0][i].size(); ++j)
            {
              std::cout << theTarOffsets[0][i][j] << std::endl;
            }
          std::cout << "src:" << i << std::endl; 
          for(unsigned j(0); j != theSrcOffsets[0][i].size(); ++j)
            {
              std::cout << theSrcOffsets[0][i][j] << std::endl;
            }
        }
      std::cout << "row1" << std::endl;
      for(unsigned i(0); i != theDiffuseSize; ++i)
        {
          std::cout << "tar:" << std::endl; 
          for(unsigned j(0); j != theTarOffsets[1][i].size(); ++j)
            {
              std::cout << theTarOffsets[1][i][j] << std::endl;
            }
          std::cout << "src:" << std::endl; 
          for(unsigned j(0); j != theSrcOffsets[1][i].size(); ++j)
            {
              std::cout << theSrcOffsets[1][i][j] << std::endl;
            }
        }
        */
    }
  void setWalkOffsets(const unsigned row, const unsigned coordA,
                      const double nDist, const Point& aLipidStart,
                      const double nLipidRadius,
                      std::vector<int>& srcOffsets)
    {
      theTarOffsets[row%2].resize(theDiffuseSize);
      theSrcOffsets[row%2].resize(theDiffuseSize);
      for(unsigned i(0); i != theDiffuseSize; ++i)
        {
          const unsigned coordB(theLattice[coordA+lipStartCoord
                               ].adjoiningCoords[i]-lipStartCoord);
          std::vector<int> tarOffsets;
          setCoordOffsets(coordB, coordA, nDist, aLipidStart, nLipidRadius,
                          tarOffsets);
          setDiffOffsets(srcOffsets, tarOffsets, theTarOffsets[row%2][i]);
          setDiffOffsets(tarOffsets, srcOffsets, theSrcOffsets[row%2][i]);
        }
    }
  void setDiffOffsets(std::vector<int>& srcOffsets,
                      std::vector<int>& tarOffsets,
                      std::vector<int>& aWalkOffsets)
    {
      for(unsigned i(0); i != tarOffsets.size(); ++i)
        {
          if(std::find(srcOffsets.begin(), srcOffsets.end(), tarOffsets[i]) ==
             srcOffsets.end())
            {
              aWalkOffsets.push_back(tarOffsets[i]);
            }
        }
    }
  void setCoordOffsets(const unsigned coordA, const unsigned coordB,
                       const double nDist, const Point& aLipidStart,
                       const double nLipidRadius,
                       std::vector<int>& anIntersectOffsets)
    {
      std::vector<unsigned> anIntersectLipids;
      getIntersectLipidsRegular(coordA+vacStartCoord, nDist, aLipidStart,
                                nLipidRadius, anIntersectLipids);
      anIntersectOffsets.resize(0);
      for(unsigned i(0); i != anIntersectLipids.size(); ++i)
        {
          anIntersectOffsets.push_back(long(anIntersectLipids[i])-long(coordB));
        }
    }
  void getIntersectLipidsRegular(const unsigned coordA, const double nDist,
                                 const Point& aLipidStart,
                                 const double nLipidRadius,
                                 std::vector<unsigned>& anIntersectLipids)
    {
      Point& pointA(*theLattice[coordA].point);
      double minY(pointA.y-aLipidStart.y-nDist*2);
      double minZ(pointA.z-aLipidStart.z-nDist*2);
      double maxY(pointA.y-aLipidStart.y+nDist*2);
      double maxZ(pointA.z-aLipidStart.z+nDist*2);
      unsigned rowStart((unsigned)std::max(minY/(nLipidRadius*sqrt(3)), 0.0));
      unsigned colStart((unsigned)std::max(minZ/(nLipidRadius*2), 0.0));
      unsigned rowEnd((unsigned)std::min(maxY/(nLipidRadius*sqrt(3)), 
                                         double(lipRows)));
      unsigned colEnd((unsigned)std::min(maxZ/(nLipidRadius*2),
                                         double(lipCols)));
      for(unsigned i(rowStart); i != rowEnd; ++i)
        {
          for(unsigned j(colStart); j != colEnd; ++j)
            {
              unsigned coord(i*lipCols+j);
              Point& pointB(*theLattice[coord+lipStartCoord].point);
              if(getDistance(&pointA, &pointB) < nDist)
                {
                  anIntersectLipids.push_back(coord);
                }
            }
        }
    }
  void setMultiscaleBindIDs(unsigned subID, unsigned prodID)
    {
      if(std::find(theMultiscaleBoundIDs.begin(), theMultiscaleBoundIDs.end(),
                   prodID) == theMultiscaleBoundIDs.end())
        {
          theMultiscaleBoundIDs.push_back(prodID);
        }
      if(std::find(theMultiscaleBindableIDs.begin(), 
                   theMultiscaleBindableIDs.end(),
                   subID) == theMultiscaleBindableIDs.end())
        {
          theMultiscaleBindableIDs.push_back(subID);
        }
      theMultiscaleBindIDs[subID] = prodID;
    }
  void setMultiscaleUnbindIDs(unsigned subID, unsigned prodID)
    {
      if(std::find(theMultiscaleBoundIDs.begin(), theMultiscaleBoundIDs.end(),
                   subID) == theMultiscaleBoundIDs.end())
        {
          theMultiscaleBoundIDs.push_back(subID);
        }
      theMultiscaleUnbindIDs[subID] = prodID;
    }
  //Get the fraction of number of nanoscopic molecules (anID) within the
  //multiscale molecule (index):
  double getMultiscaleBoundFraction(unsigned index, unsigned anID)
    {
      double fraction(0);
      if(isMultiscale)
        {
          const unsigned coordA(getCoord(index)-vacStartCoord);
          //Need to optimize this, which is a big bottleneck:
          if(isRegularLattice)
            {
              const int rowA(coordA/lipCols);
              const std::vector<int>& anOffsetsA(theOffsets[rowA%2]);
              unsigned size(0);
              for(unsigned i(0); i != anOffsetsA.size(); ++i)
                {
                  const int offsetRow((anOffsetsA[i]+theRegLatticeCoord)/
                                      lipCols-theRegLatticeCoord/lipCols);
                  int coord(coordA+anOffsetsA[i]);
                  if(isInLattice(coord, offsetRow+rowA))
                    {
                      ++size;
                      if(theLattice[coord+lipStartCoord].id == anID)
                        {
                          fraction += 1;
                        }
                    }
                }
              fraction /= size;
            }
          else
            {
              for(unsigned i(0); i != theIntersectLipids[coordA].size(); ++i)
                {
                  unsigned aCoord(theIntersectLipids[coordA][i]+lipStartCoord);
                  if(theLattice[aCoord].id == anID)
                    {
                      fraction += 1;
                    }
                }
              fraction /= theIntersectLipids[coordA].size();
            }
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
              const unsigned coordA(getCoord(i)-vacStartCoord);
              if(isRegularLattice)
                {
                  const int rowA(coordA/lipCols);
                  const std::vector<int>& anOffsetsA(theOffsets[
                                                     rowA%2]);
                  for(unsigned j(0); j != anOffsetsA.size(); ++j)
                    {
                      const int offsetRow((anOffsetsA[j]+theRegLatticeCoord)/
                                          lipCols-theRegLatticeCoord/lipCols);
                      int coord(coordA+anOffsetsA[j]);
                      if(isInLattice(coord, offsetRow+rowA))
                        {
                          coord += lipStartCoord;
                          if(theLattice[coord].id == theID)
                            {
                              thePopulatableCoords.push_back(coord);
                            }
                        }
                    }
                }
              else
                {
                  for(unsigned j(0); j != theIntersectLipids[coordA].size();
                      ++j)
                    {
                      unsigned aCoord(theIntersectLipids[coordA][j]+
                                      lipStartCoord);
                      if(theLattice[aCoord].id == theID)
                        {
                          thePopulatableCoords.push_back(aCoord);
                        }
                    }
                }
            }
          return thePopulatableCoords.size();
        }
      //Required by populate dense because some comp vacant voxels would
      //have become interface species voxels and no longer populatable:
      else if(isCompVacant)
        {
          unsigned aSize(0);
          for(unsigned i(0); i != theMoleculeSize; ++i)
            {
              if(theMolecules[i]->id == theID)
                {
                  ++aSize;
                }
            }
          return aSize;
        }
      return theMoleculeSize;
    }
  Voxel* getRandomPopulatableMolecule()
    {
      Voxel* aMolecule;
      if(isMultiscale)
        {
          unsigned index(0);
          do
            {
              index = gsl_rng_uniform_int(theRng, thePopulatableCoords.size());
            }
          while(theLattice[thePopulatableCoords[index]].id != theID);
          aMolecule =  &theLattice[thePopulatableCoords[index]];
        }
      else
        {
          aMolecule = getRandomMolecule();
          while(aMolecule->id != theID)
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
  bool isPopulatable(Voxel* aVoxel)
    {
      if(isMultiscale)
        {
          if(isIntersectMultiscale(aVoxel->coord))
            {
              return false;
            }
        }
      else if(aVoxel->id != theVacantID)
        {
          return false;
        }
      return true;
    }
  //Can aVoxel of this species replaced by aSpecies:
  bool isReplaceable(Voxel* aVoxel, Species* aSpecies)
    {
      if(getComp() != aSpecies->getComp() &&
         theID != aSpecies->getVacantID())
        {
          return false;
        }
      if(aSpecies->getIsMultiscale())
        {
          if(aSpecies->isIntersectMultiscale(aVoxel->coord))
            {
              return false;
            }
        }
      if(getVacantID() != aSpecies->getID() && 
         aSpecies->getVacantID() != theID && 
         getVacantID() != aSpecies->getVacantID())
        {
          return false;
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
  bool isPeriodic;
  bool isPolymer;
  bool isReactiveVacant;
  bool isRegularLattice;
  bool isSubunitInitialized;
  bool isTag;
  bool isTagged;
  bool isVacant;
  const unsigned short theID;
  int lipCols;
  int lipRows;
  unsigned lipStartCoord;
  unsigned theAdjoiningCoordSize;
  unsigned theCollision;
  unsigned theDiffuseSize;
  unsigned theDimension;
  unsigned theInitCoordSize;
  unsigned theMoleculeSize;
  unsigned theNullCoord;
  unsigned theNullID;
  unsigned theRegLatticeCoord;
  unsigned theSpeciesSize;
  unsigned vacCols;
  unsigned vacRows;
  unsigned vacStartCoord;
  int thePolymerDirectionality;
  int theVacantID;
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
  std::vector<std::vector<int> > theOffsets;
  std::vector<std::vector<std::vector<int> > > theTarOffsets;
  std::vector<std::vector<std::vector<int> > > theSrcOffsets;
  std::vector<bool> theFinalizeReactions;
  std::vector<unsigned> collisionCnts;
  std::vector<unsigned> theCoords;
  std::vector<unsigned> theMultiscaleBindIDs;
  std::vector<unsigned> theMultiscaleBoundIDs;
  std::vector<unsigned> theMultiscaleUnbindIDs;
  std::vector<unsigned> theMultiscaleBindableIDs;
  std::vector<unsigned> thePopulatableCoords;
  std::vector<Tag> theTags;
  std::vector<double> theBendAngles;
  std::vector<double> theReactionProbabilities;
  std::vector<Voxel*> theMolecules;
  std::vector<Voxel*>* theCompVoxels;
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

