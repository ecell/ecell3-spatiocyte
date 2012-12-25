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
#include <algorithm>
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
 * that theMols list is not updated when other species diffuse on it. 
 * There are four possible types of vacant species:
 * 1. !isDiffusiveVacant && !isReactiveVacant: a vacant species
 *    whose theMols list and theMolSize are never updated. This is
 *    the most basic version. This vacant species never diffuses, or reacts
 *    using SNRP. 
 * 2. !isDiffusiveVacant && isReactiveVacant: a vacant species which is a
 *    substrate of a SNRP reaction. In this case theMolSize is updated
 *    before the step interval of SNRP is calculated, and theMols list is
 *    updated before the SNRP reaction (fire) is executed to get a valid
 *    molecule list. Set by SNRP.
 * 3. isDiffusiveVacant && !isReactiveVacant: a vacant species which also
 *    diffuses. In this case, theMols list and theMolSize are
 *    updated just before it is diffused. Set by DiffusionProcess.   
 * 4. isDiffusiveVacant && isReactiveVacant: a vacant species which reacts
 *    using SNRP and also diffuses.
 * 5. isCompVacant: the VACANT species declared for each compartment. It can
 *    also be isDiffusiveVacant and isReactiveVacant. Set during compartment
 *    registration. It also persistently stores all the compartment voxels.
 *    Referred to as theVacantSpecies. For isCompVacant, initially there
 *    are no molecules in its list. All voxels are stored in theCompVoxels. Only
 *    if it is updated before being called by VisualizationLogProcess, SNRP or
 *    DiffusionProcess, it will be theMols will be populated with the
 *    CompMols.
 * 6. isVacant {isCompVacant; isDiffusiveVacant; isReactiveVacant): the
 *    general name used to identify either isCompVacant, isDiffusiveVacant or
 *    isReactiveVacant to reduce comparison operations.
 */

class Species
{
public:
  Species(SpatiocyteStepper* aStepper, Variable* aVariable,
          unsigned anID, int anInitMolSize, const gsl_rng* aRng,
          double voxelRadius, std::vector<unsigned>& aLattice,
          std::vector<VoxelInfo>& anInfo):
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
    theInitMolSize(anInitMolSize),
    theMolSize(0),
    D(0),
    theDiffuseRadius(voxelRadius),
    theDiffusionInterval(libecs::INF),
    theMolRadius(voxelRadius),
    theVoxelRadius(voxelRadius),
    theWalkProbability(1),
    theRng(aRng),
    thePopulateProcess(NULL),
    theStepper(aStepper),
    theVariable(aVariable),
    theCompMols(&theMols),
    theInfo(anInfo),
    theLattice(aLattice) {}
  ~Species() {}
  void initialize(int speciesSize, int anAdjoinSize,
                  unsigned aNullMol, unsigned aNullID)
    {
      theAdjoinSize = anAdjoinSize;
      theNullMol = aNullMol;
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
      theNullTag.origin = theNullMol;
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
      else if(theMolSize)
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
      else if(theMolSize)
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
      else if(theMolSize)
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
      else if(theMolSize)
        {
          std::cout << "Species:" << theVariable->getFullID().asString() <<
            " not MoleculePopulated." << std::endl;
        }
    }
  void populateUniformOnMultiscale()
    {
      if(thePopulateProcess)
        {
          thePopulateProcess->populateUniformOnMultiscale(this);
        }
      else if(theMolSize)
        {
          std::cout << "Species:" << theVariable->getFullID().asString() <<
            " not MoleculePopulated." << std::endl;
        }
    }
  Variable* getVariable() const
    {
      return theVariable;
    }
  unsigned size() const
    {
      return theMolSize;
    }
  unsigned& getMol(unsigned anIndex)
    {
      return theMols[anIndex];
    }
  Point& getPoint(unsigned anIndex)
    {
      return theInfo[getMol(anIndex)].point;
    }
  unsigned getID() const
    {
      return theID;
    }
  double getMeanSquaredDisplacement()
    {
      if(!theMolSize)
        {
          return 0;
        }
      double aDisplacement(0);
      for(unsigned i(0); i < theMolSize; ++i)
        {
          Point aCurrentPoint(theStepper->getPeriodicPoint(
                                                 getMol(i),
                                                 theDimension,
                                                 &theMolOrigins[i]));
          double aDistance(getDistance(theMolOrigins[i].point,
                                       aCurrentPoint));
          aDisplacement += aDistance*aDistance;
        }
      return aDisplacement*pow(theDiffuseRadius*2, 2)/theMolSize;
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
      theInitMolSize = theMolSize;
      getVariable()->setValue(theMolSize);
    }
  void finalizeSpecies()
    {
      if(theCollision)
        {
          collisionCnts.resize(theMolSize);
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
                  std::random_shuffle(theMols.begin(), theMols.end());
                  break;
                }
            }
          if(isTagged)
            {
              for(unsigned i(0); i != theMolSize; ++i)
                {
                  theTags[i].origin = getMol(i);
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
      return theMolSize == theInitMolSize;
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
  void addCollision(unsigned aMol)
    {
      for(unsigned i(0); i < theMolSize; ++i)
        {
          if(theMols[i] == aMol)
            {
              ++collisionCnts[i];
              return;
            }
        }
      std::cout << "error in species add collision" << std::endl;
    }
  void setTarMols()
    {
      if(theTarMols.size() < theMolSize)
        {
          theTarMols.resize(theMolSize);
          theRands.resize(std::min(unsigned(10000), theMolSize*10));
          for(unsigned i(0); i != theRands.size(); ++i)
            {
              theRands[i] = gsl_rng_uniform_int(theRng, theAdjoinSize-1)+1;
            }
        }
      /*
      unsigned j(gsl_rng_uniform_int(theRng, theRands.size()));
      for(unsigned i(0); i != theMolSize; ++i, ++j)
        {
          while(j == theRands.size())
            {
              j = gsl_rng_uniform_int(theRng, theRands.size());
            }
          theTarMols[i] = theLattice[theMols[i]*theAdjoinSize+theRands[j]];
        }
        */
    }
  void walk()
    {
      const unsigned beginMolSize(theMolSize);
      setTarMols();
      unsigned j(gsl_rng_uniform_int(theRng, theRands.size()));
      for(unsigned i(0); i < beginMolSize && i < theMolSize; ++i, ++j)
        {
          while(j == theRands.size())
            {
              j = gsl_rng_uniform_int(theRng, theRands.size());
            }
          const unsigned tar();
          if(theLattice[theLattice[theMols[i]*theAdjoinSize+theRands[j]]*theAdjoinSize] == theVacantID)
            {
              if(theWalkProbability == 1 ||
                 gsl_rng_uniform(theRng) < theWalkProbability)
                {
                  theLattice[theLattice[theMols[i]*theAdjoinSize+theRands[j]]*theAdjoinSize] = theID;
                  theLattice[theMols[i]*theAdjoinSize] = theVacantID;
                  theMols[i] = theLattice[theMols[i]*theAdjoinSize+theRands[j]];
                }
            }
          /*
          else
            {
              unsigned targetMol(theTarMols[i]);
              if(theLattice[theTarMols[i]].id == theComp->interfaceID)
                {
                  unsigned diffuseSize(theLattice[targetMol].diffuseSize);
                  unsigned range(theInfo[targetMol].adjoinSize-diffuseSize);
                  unsigned index(gsl_rng_uniform_int(theRng, range));
                  targetMol = theLattice[targetMol*theAdjoinSize+diffuseSize+index];
                }
              const unsigned& targetID(theLattice[targetMol].id);
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
                          targetIndex = targetSpecies->getMolIndex(targetMol);
                        }
                      if(theCollision)
                        { 
                          ++collisionCnts[i];
                          targetSpecies->addCollision(targetMol);
                          if(theCollision != 2)
                            {
                              return;
                            }
                        }
                      unsigned aMolSize(theMolSize);
                      react(theMols[i], targetMol, i, targetIndex,
                            targetSpecies);
                      //If the reaction is successful, the last molecule of this
                      //species will replace the pointer of i, so we need to 
                      //decrement i to perform the diffusion on it. However, if
                      //theMolSize didn't decrease, that means the
                      //currently walked molecule was a product of this
                      //reaction and so we don't need to walk it again by
                      //decrementing i.
                      if(theMolSize < aMolSize)
                        {
                          --i;
                        }
                    }
                }
            }
            */
        }
    }
  void walkMultiscale()
    {
      unsigned beginMolSize(theMolSize);
      for(unsigned i(0); i < beginMolSize && i < theMolSize; ++i)
        {
          unsigned srcMol(theMols[i]);
          unsigned& source(theLattice[srcMol*theAdjoinSize]);
          int size(theInfo[srcMol].diffuseSize);
          unsigned tarMol(theLattice[theMols[i]*theAdjoinSize+gsl_rng_uniform_int(theRng, size-1)+1]);
          unsigned& target(theLattice[tarMol*theAdjoinSize]);
          if(target == theVacantID)
            {
              if(theWalkProbability == 1 ||
                 gsl_rng_uniform(theRng) < theWalkProbability)
                {
                  if(!isIntersectMultiscale(srcMol, tarMol))
                    {
                      removeMultiscaleMol(srcMol);
                      addMultiscaleMol(tarMol);
                      target = theID;
                      source = theVacantID;
                      theMols[i] = tarMol;
                    }
                }
            }
        }
    }
  void walkVacant()
    {
      updateVacantMols();
      for(unsigned i(0); i < theMolSize; ++i)
        {
          unsigned& source(theLattice[theMols[i]*theAdjoinSize]);
          int size;
          if(isFixedAdjoins)
            {
              size = theAdjoinSize;
            }
          else
            {
              size = theInfo[theMols[i]].diffuseSize;
            }
          unsigned aMol(theLattice[theMols[i]*theAdjoinSize+
                        gsl_rng_uniform_int(theRng, size-1)+1]);
          unsigned& target(theLattice[aMol*theAdjoinSize]);
          if(target == theVacantID)
            {
              if(theWalkProbability == 1 ||
                 gsl_rng_uniform(theRng) < theWalkProbability)
                {
                  target = theID;
                  source = theVacantID;
                  theMols[i] = aMol;
                }
            }
        }
    }
  void react(unsigned srcMol, unsigned tarMol, unsigned sourceIndex,
             unsigned targetIndex, Species* targetSpecies)
    {
      DiffusionInfluencedReactionProcess* aReaction(
               theDiffusionInfluencedReactions[targetSpecies->getID()]);
      unsigned moleculeA(srcMol);
      unsigned moleculeB(tarMol);
      unsigned indexA(sourceIndex);
      unsigned indexB(targetIndex);
      if(aReaction->getA() != this)
        {
          indexA = targetIndex; 
          indexB = sourceIndex;
          moleculeA = tarMol;
          moleculeB = srcMol;
        }
      if(aReaction->react(moleculeA, moleculeB, indexA, indexB))
        {
          //Soft remove the source molecule, i.e., keep the id:
          softRemoveMolIndex(sourceIndex);
          //Soft remove the target molecule:
          //Make sure the targetIndex is valid:
          //Target and Source are same species:
          //For some reason if I use theMols[sourceIndex] instead
          //of getMol(sourceIndex) the walk method becomes
          //much slower when it is only diffusing without reacting:
          if(targetSpecies == this && getMol(sourceIndex) == tarMol)
            {
              softRemoveMolIndex(sourceIndex);
            }
          //If the targetSpecies is a multiscale species with implicit
          //molecule, theTargetIndex is equal to the target molecule size,
          //so we use this info to avoid removing the implicit target molecule:
          else if(targetIndex != targetSpecies->size())
            {
              targetSpecies->softRemoveMolIndex(targetIndex);
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
      int newMolSize(0);
      for(unsigned i(0); i < theMolSize; ++i) 
        {
          unsigned aMol(getMol(i));
          if(theStepper->isRemovableEdgeMol(aMol, theComp))
            {
              Comp* aSuperComp(
                 theStepper->system2Comp(theComp->system->getSuperSystem())); 
              aSuperComp->vacantSpecies->addCompVoxel(aMol);
            }
          else 
            { 
              theMols[newMolSize] = aMol;
              ++newMolSize; 
            }
        }
      theMolSize = newMolSize;
      //Must resize, otherwise compVoxelSize will be inaccurate:
      theMols.resize(theMolSize);
      theVariable->setValue(theMolSize);
    }
  void removePeriodicEdgeVoxels()
    {
      int newMolSize(0);
      for(unsigned i(0); i < theMolSize; ++i) 
        {
          unsigned aMol(getMol(i));
          if(theStepper->isPeriodicEdgeMol(aMol, theComp))
            {
              theLattice[aMol*theAdjoinSize] = theLattice[theNullMol*theAdjoinSize];
            }
          else 
            { 
              theMols[newMolSize] = aMol;
              ++newMolSize; 
            }
        }
      theMolSize = newMolSize;
      //Must resize, otherwise compVoxelSize will be inaccurate:
      theMols.resize(theMolSize);
      theVariable->setValue(theMolSize);
    }
  void updateSpecies()
    {
      if(isCompVacant && (isDiffusiveVacant || isReactiveVacant))
        {
          theCompMols = new std::vector<unsigned>;
          for(unsigned i(0); i != theMolSize; ++i)
            { 
              theCompMols->push_back(theMols[i]);
            }
        }
    }
  //If it isReactiveVacant it will only be called by SNRP when it is substrate
  //If it isDiffusiveVacant it will only be called by DiffusionProcess before
  //being diffused. So we need to only check if it isVacant:
  void updateMols()
    {
      if(isDiffusiveVacant || isReactiveVacant)
        {
          updateVacantMols();
        }
      else if(isTag)
        {
          updateTagMols();
        }
    }
  //If it isReactiveVacant it will only be called by SNRP when it is substrate:
  void updateMolSize()
    {
      if(isDiffusiveVacant || isReactiveVacant)
        {
          updateVacantMolSize();
        }
    }
  void updateTagMols()
    {
      theMolSize = 0;
      for(unsigned i(0); i != theTaggedSpeciesList.size(); ++i)
        {
          Species* aSpecies(theTaggedSpeciesList[i]);
          for(unsigned j(0); j != aSpecies->size(); ++j)
            {
              if(aSpecies->getTagID(j) == theID)
                {
                  unsigned aMol(aSpecies->getMol(j));
                  ++theMolSize;
                  if(theMolSize > theMols.size())
                    {
                      theMols.push_back(aMol);
                    }
                  else
                    {
                      theMols[theMolSize-1] = aMol;
                    }
                }
            }
        }
    }
  //Even if it is a isCompVacant, this method will be called by
  //VisualizationLogProcess, or SNRP if it is Reactive, or DiffusionProcess
  //if it is Diffusive:
  void updateVacantMols()
    {
      theMolSize = 0;
      int aSize(theVacantSpecies->compVoxelSize());
      for(int i(0); i != aSize; ++i)
        { 
          //Voxel* aVoxel(theVacantSpecies->getCompVoxel(i));
          unsigned aMol(theVacantSpecies->getCompMol(i));
          //if(aVoxel->id == theID)
          if(theLattice[aMol*theAdjoinSize] == theID)
            {
              ++theMolSize;
              if(theMolSize > theMols.size())
                {
                  theMols.push_back(aMol);
                }
              else
                {
                  theMols[theMolSize-1] = aMol;
                }
            }
        }
      theVariable->setValue(theMolSize);
    }
  void updateVacantMolSize()
    {
      theMolSize = 0;
      int aSize(theVacantSpecies->compVoxelSize());
      for(int i(0); i != aSize; ++i)
        { 
          //Voxel* aVoxel(theVacantSpecies->getCompVoxel(i));
          unsigned aMol(theVacantSpecies->getCompMol(i));
          //if(aVoxel->id == theID)
          if(theLattice[aMol*theAdjoinSize] == theID)
            {
              ++theMolSize;
            }
        }
      if(theMolSize > theMols.size())
        {
          theMols.resize(theMolSize);
        }
      theVariable->setValue(theMolSize);
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
      if(isTagged && anIndex != theMolSize)
        {
          return theTags[anIndex];
        }
      return theNullTag;
    }
  void addMol(unsigned aMol)
    {
      addMol(aMol, theNullTag);
    }
  Species* getMultiscaleVacantSpecies()
    {
      return theMultiscaleVacantSpecies;
    }
  void addMol(unsigned aMol, Tag& aTag)
    {
      if(isMultiscale)
        {
          Species* aSpecies(theStepper->id2species(theLattice[aMol*theAdjoinSize]));
          if(aSpecies->getVacantSpecies() != theMultiscaleVacantSpecies)
            {
              doAddMol(aMol, aTag);
              addMultiscaleMol(aMol);
            }
        }
      else if(!isVacant)
        {
          doAddMol(aMol, aTag);
        }
      theLattice[aMol*theAdjoinSize] = theID;
    }
  void doAddMol(unsigned aMol, Tag& aTag)
    {
      ++theMolSize; 
      if(theMolSize > theMols.size())
        {
          theMols.resize(theMolSize);
          theTags.resize(theMolSize);
        }
      theMols[theMolSize-1] = aMol;
      if(isTagged)
        {
          //If it is theNullTag:
          if(aTag.origin == theNullMol)
            {
              Tag aNewTag = {getMol(theMolSize-1), theNullID};
              theTags[theMolSize-1] = aNewTag;
            }
          else
            {
              theTags[theMolSize-1] = aTag;
            }
        }
      theVariable->setValue(theMolSize);
    }
  void addMultiscaleMol(unsigned aMol)
    {
      unsigned coordA(aMol-vacStartMol);
      for(unsigned i(0); i != theIntersectLipids[coordA].size(); ++i)
        {
          unsigned coordB(theIntersectLipids[coordA][i]+lipStartMol);
          if(theLattice[coordB*theAdjoinSize] == theMultiscaleVacantSpecies->getID())
            {
              theLattice[coordB*theAdjoinSize] = theID;
            }
          else
            {
              multiscaleBind(coordB);
            }
        }
    }
  void removeMultiscaleMol(unsigned aMol)
    {
      unsigned coordA(aMol-vacStartMol);
      for(unsigned i(0); i != theIntersectLipids[coordA].size(); ++i)
        {
          unsigned coordB(theIntersectLipids[coordA][i]+lipStartMol);
          if(theLattice[coordB*theAdjoinSize] == theID)
            {
              theLattice[coordB*theAdjoinSize] = theMultiscaleVacantSpecies->getID();
            }
          else
            {
              multiscaleUnbind(coordB);
            }
        }
    }
  bool isIntersectMultiscale(unsigned aMol)
    {
      unsigned coordA(aMol-vacStartMol);
      for(unsigned i(0); i != theIntersectLipids[coordA].size(); ++i)
        {
          unsigned coordB(theIntersectLipids[coordA][i]+lipStartMol);
          unsigned anID(theLattice[coordB*theAdjoinSize]);
          if(anID == theID ||
             std::find(theMultiscaleIDs.begin(), theMultiscaleIDs.end(),
                       anID) != theMultiscaleIDs.end())
            {
              return true;
            }
        }
      return false;
    }
  bool isIntersectMultiscale(unsigned srcMol, unsigned tarMol)
    {
      bool isIntersect(false);
      unsigned coordA(srcMol-vacStartMol);
      std::vector<unsigned> temp;
      temp.resize(theIntersectLipids[coordA].size());
      for(unsigned i(0); i != theIntersectLipids[coordA].size(); ++i)
        {
          unsigned coordB(theIntersectLipids[coordA][i]+lipStartMol);
          temp[i] = theLattice[coordB*theAdjoinSize];
          theLattice[coordB*theAdjoinSize] = theSpeciesSize;
        }
      isIntersect = isIntersectMultiscale(tarMol);
      coordA = srcMol-vacStartMol;
      for(unsigned i(0); i != theIntersectLipids[coordA].size(); ++i)
        {
          unsigned coordB(theIntersectLipids[coordA][i]+lipStartMol);
          theLattice[coordB*theAdjoinSize] = temp[i];
        }
      return isIntersect;
    }
  void multiscaleBind(unsigned aMol)
    {
      unsigned anID(theLattice[aMol*theAdjoinSize]);
      Species* source(theStepper->id2species(anID));
      Species* target(theStepper->id2species(theMultiscaleBindIDs[anID]));
      source->softRemoveMol(aMol);
      target->addMol(aMol);
    }
  void multiscaleUnbind(unsigned aMol)
    {
      unsigned anID(theLattice[aMol*theAdjoinSize]);
      Species* source(theStepper->id2species(anID));
      Species* target(theStepper->id2species(theMultiscaleUnbindIDs[anID]));
      source->softRemoveMol(aMol);
      target->addMol(aMol);
    }
  void addCompVoxel(unsigned aMol)
    {
      theLattice[aMol*theAdjoinSize] = theID;
      //theCompVoxels->push_back(&theLattice[aMol]);
      theCompMols->push_back(aMol);
      ++theMolSize;
      theVariable->setValue(theMolSize);
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
      return theCompMols->size();
    }
  unsigned getCompMol(unsigned index)
    {
      return (*theCompMols)[index];
    }
  /*
  Voxel* getCompVoxel(unsigned index)
    {
      return (*theCompVoxels)[index];
    }
    */
  unsigned getMolIndexFast(unsigned aMol)
    {
      for(unsigned i(0); i < theMolSize; ++i)
        {
          if(theMols[i] == aMol)
            {
              return i;
            }
        }
      return theMolSize;
    }
  unsigned getMolIndex(unsigned aMol)
    {
      unsigned index(getMolIndexFast(aMol));
      if(index == theMolSize)
        { 
          if(isDiffusiveVacant || isReactiveVacant)
            {
              updateVacantMols();
            }
          index = getMolIndexFast(aMol);
          if(index == theMolSize)
            { 
              std::cout << "error in getting the index:" << getIDString() <<
               " size:" << theMolSize << std::endl;
              return 0;
            }
        }
      return index;
    }
  //it is soft remove because the id of the molecule is not changed:
  void softRemoveMol(unsigned aMol)
    {
      if(isMultiscale)
        {
          Species* aSpecies(theStepper->id2species(theLattice[aMol*theAdjoinSize]));
          if(aSpecies->getVacantSpecies() != theMultiscaleVacantSpecies)
            {
              softRemoveMolIndex(getMolIndex(aMol));
            }
        }
      else if(!isVacant)
        {
          softRemoveMolIndex(getMolIndex(aMol));
        }
    }
  void removeMol(unsigned aMol)
    {
      if(isMultiscale)
        {
          Species* aSpecies(theStepper->id2species(theLattice[aMol*theAdjoinSize]));
          if(aSpecies->getVacantSpecies() != theMultiscaleVacantSpecies)
            {
              removeMolIndex(getMolIndex(aMol));
            }
        }
      if(!isVacant)
        {
          removeMolIndex(getMolIndex(aMol));
        }
    }
  void removeMolIndex(unsigned anIndex)
    {
      if(!isVacant)
        {
          theLattice[theMols[anIndex]*theAdjoinSize] = theVacantID;
          softRemoveMolIndex(anIndex);
        }
    }
  void softRemoveMolIndex(unsigned anIndex)
    {
      if(isMultiscale)
        {
          removeMultiscaleMol(theMols[anIndex]);
        }
      if(!isVacant)
        {
          theMols[anIndex] = theMols[--theMolSize];
          if(isTagged)
            {
              theTags[anIndex] = theTags[theMolSize];
            }
          theVariable->setValue(theMolSize);
          return;
        }
    }
  //Used to remove all molecules and free memory used to store the molecules
  void clearMols()
    {
      theMols.resize(0);
      theMolSize = 0;
      theVariable->setValue(0);
    }
  //Used by the SpatiocyteStepper when resetting an interation, so must
  //clear the whole compartment using theComp->vacantSpecies->getVacantID():
  void removeMols()
    {
      if(!isCompVacant)
        {
          for(unsigned i(0); i < theMolSize; ++i)
            {
              theLattice[theMols[i]*theAdjoinSize] = theVacantSpecies->getID();
            }
          theMolSize = 0;
          theVariable->setValue(theMolSize);
        }
    }
  int getPopulateMolSize()
    {
      return theInitMolSize-theMolSize;
    }
  int getInitMolSize()
    {
      return theInitMolSize;
    }
  void initMolOrigins()
    {
      theMolOrigins.resize(theMolSize);
      for(unsigned i(0); i < theMolSize; ++i)
        {
          Origin& anOrigin(theMolOrigins[i]);
          anOrigin.point = theStepper->coord2point(getMol(i));
          anOrigin.row = 0;
          anOrigin.layer = 0;
          anOrigin.col = 0;
        }
    }
  void removeBoundaryMols()
    {
      for(unsigned i(0); i < theMolSize; ++i)
        {
          if(theStepper->isBoundaryMol(getMol(i), theDimension))
            {
              std::cout << "is still there" << std::endl;
            }
        }
      theVariable->setValue(theMolSize);
    }
  void relocateBoundaryMols()
    {
      for(unsigned i(0); i < theMolSize; ++i)
        {
          Origin anOrigin(theMolOrigins[i]);
          unsigned periodicMol(theStepper->getPeriodicMol(
                                                getMol(i),
                                                theDimension, &anOrigin));
          if(theLattice[periodicMol*theAdjoinSize] == theVacantID)
            {
              theLattice[theMols[i]*theAdjoinSize] = theVacantID;
              theMols[i] = periodicMol;
              theLattice[theMols[i]*theAdjoinSize] = theID;
              theMolOrigins[i] = anOrigin;
            }
        }
    }
  unsigned getVacantID() const
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
  double getMolRadius() const
    {
      return theMolRadius;
    }
  double getDiffuseRadius() const
    {
      return theDiffuseRadius;
    }
  void setMolRadius(double aRadius)
    {
      theMolRadius = aRadius;
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
      return gsl_rng_uniform_int(theRng, theMolSize);
    }
  unsigned getRandomMol()
    {
      return theMols[getRandomIndex()];
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
  unsigned getRandomAdjoin(unsigned srcMol, int searchVacant)
    {
      std::vector<unsigned> compMols;
      if(searchVacant)
        { 
          for(unsigned i(1); i != theInfo[srcMol].adjoinSize; ++i)
            {
              unsigned aMol(theLattice[srcMol*theAdjoinSize+i]);
              if(isPopulatable(aMol))
                {
                  compMols.push_back(aMol);
                }
            }
        }
      else
        {
          for(unsigned i(1); i != theInfo[srcMol].adjoinSize; ++i)
            {
              unsigned aMol(theLattice[srcMol*theAdjoinSize+i]);
              if(theStepper->id2Comp(theLattice[aMol*theAdjoinSize]) == theComp)
                {
                  compMols.push_back(aMol);
                }
            }
        }
      return getRandomVacantMol(compMols);
    } 
  unsigned getBindingSiteAdjoin(unsigned srcMol, int bindingSite)
    {
      if(bindingSite < theInfo[srcMol].adjoinSize)
        { 
          unsigned aMol(theLattice[srcMol*theAdjoinSize+bindingSite+1]);
          if(isPopulatable(aMol))
            {
              return aMol;
            }
        }
      return theNullMol;
    } 
  unsigned getRandomAdjoin(unsigned srcMol, Species* aTargetSpecies,
                           int searchVacant)
    {
      std::vector<unsigned> compMols;
      if(searchVacant)
        { 
          for(unsigned i(1); i != theInfo[srcMol].adjoinSize; ++i)
            {
              unsigned aMol(theLattice[srcMol*theAdjoinSize+i]);
              if(theLattice[aMol*theAdjoinSize] == aTargetSpecies->getID())
                {
                  compMols.push_back(aMol);
                }
            }
        }
      else
        {
          for(unsigned i(1); i != theInfo[srcMol].adjoinSize; ++i)
            {
              unsigned aMol(theLattice[srcMol*theAdjoinSize+i]);
              if(theStepper->id2Comp(theLattice[aMol*theAdjoinSize]) == theComp)
                {
                  compMols.push_back(aMol);
                }
            }
        }
      return getRandomVacantMol(compMols, aTargetSpecies);
    } 
  unsigned getRandomAdjoin(unsigned srcMol, unsigned tarMol,
                           int searchVacant)
    {
      std::vector<unsigned> compMols;
      if(searchVacant)
        { 
          for(unsigned i(1); i != theInfo[srcMol].adjoinSize; ++i)
            {
              unsigned aMol(theLattice[srcMol*theAdjoinSize+i]);
              if(aMol != tarMol && isPopulatable(aMol))
                {
                  compMols.push_back(aMol);
                }
            }
        }
      else
        {
          for(unsigned i(1); i != theInfo[srcMol].adjoinSize; ++i)
            {
              unsigned aMol(theLattice[srcMol*theAdjoinSize+i]);
              if(theStepper->id2Comp(theLattice[aMol*theAdjoinSize]) == theComp &&
                 aMol != tarMol)
                {
                  compMols.push_back(aMol);
                }
            }
        }
      return getRandomVacantMol(compMols);
    }
  unsigned getAdjoinMolCnt(unsigned srcMol, Species* aTargetSpecies)
    {
      unsigned cnt(0);
      for(unsigned i(1); i != theInfo[srcMol].adjoinSize; ++i)
        {
          if(theLattice[theLattice[srcMol*theAdjoinSize+i]*theAdjoinSize] == aTargetSpecies->getID())
            {
              ++cnt;
            }
        }
      return cnt;
    }
  unsigned getRandomAdjoin(unsigned srcMol, unsigned targetA,
                           unsigned targetB, int searchVacant)
    {
      std::vector<unsigned> compMols;
      if(srcMol != targetA && srcMol != targetB)
        {
          if(searchVacant)
            { 
              for(unsigned i(1); i != theInfo[srcMol].adjoinSize; ++i)
                {
                  unsigned aMol(theLattice[srcMol*theAdjoinSize+i]);
                  if(isPopulatable(aMol))
                    {
                      compMols.push_back(aMol);
                    }
                }
            }
          else
            {
              for(unsigned i(1); i != theInfo[srcMol].adjoinSize; ++i)
                {
                  unsigned aMol(theLattice[srcMol*theAdjoinSize+i]);
                  if(theStepper->id2Comp(theLattice[aMol]*theAdjoinSize) == theComp)
                    {
                      compMols.push_back(aMol);
                    }
                }
            }
        }
      return getRandomVacantMol(compMols);
    }
  unsigned getRandomVacantMol(std::vector<unsigned>& aMols)
    {
      if(aMols.size())
        {
          const int r(gsl_rng_uniform_int(theRng, aMols.size())); 
          unsigned aMol(aMols[r]);
          if(isPopulatable(aMol))
            {
              return aMol;
            }
        }
      return theNullMol;
    }
  unsigned getRandomVacantMol(std::vector<unsigned>& aMols,
                                Species* aVacantSpecies)
    {
      if(aMols.size())
        {
          const int r(gsl_rng_uniform_int(theRng, aMols.size())); 
          unsigned aMol(aMols[r]);
          if(theLattice[aMol*theAdjoinSize] == aVacantSpecies->getID())
            {
              return aMol;
            }
        }
      return theNullMol;
    }
  unsigned getRandomCompMol(int searchVacant)
    {
      Species* aVacantSpecies(theComp->vacantSpecies);
      int aSize(aVacantSpecies->compVoxelSize());
      int r(gsl_rng_uniform_int(theRng, aSize));
      if(searchVacant)
        {
          for(int i(r); i != aSize; ++i)
            {
              unsigned aMol(aVacantSpecies->getCompMol(i));
              if(isPopulatable(aMol))
                {
                  return aMol;
                }
            }
          for(int i(0); i != r; ++i)
            {
              unsigned aMol(aVacantSpecies->getCompMol(i));
              if(isPopulatable(aMol))
                {
                  return aMol;
                }
            }
        }
      else
        {
          unsigned aMol(aVacantSpecies->getCompMol(r));
          if(isPopulatable(aMol))
            {
              return aMol;
            }
        }
      return theNullMol;
    }
  unsigned getRandomAdjoinCompMol(Comp* aComp, int searchVacant)
    {
      int aSize(theVacantSpecies->size());
      int r(gsl_rng_uniform_int(theRng, aSize)); 
      return getRandomAdjoin(theVacantSpecies->getMol(r), searchVacant);
    }
  void setVacStartMol(unsigned aMol)
    {
      vacStartMol = aMol;
    }
  void setLipStartMol(unsigned aMol)
    {
      lipStartMol = aMol;
    }
  void setIntersectLipids(Species* aLipid)
    {
      //Traverse through the entire compartment voxels:
      unsigned endA(vacStartMol+theVacantSpecies->size());
      unsigned endB(lipStartMol+aLipid->size());
      double dist((aLipid->getMolRadius()+theMolRadius)/
                  (2*theVoxelRadius));
      theIntersectLipids.resize(theVacantSpecies->size());
      for(unsigned i(vacStartMol); i != endA; ++i)
        {
          Point& pointA(theInfo[i].point);
          for(unsigned j(lipStartMol); j != endB; ++j)
            {
              Point& pointB(theInfo[j].point);
              if(getDistance(pointA, pointB) < dist)
                {
                  //We save j-lipStartMol and not the absolute coord
                  //since the absolute coord may change after resize 
                  //of lattice:
                  theIntersectLipids[i-vacStartMol
                    ].push_back(j-lipStartMol);
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
          unsigned i(getMol(index)-vacStartMol);
          for(unsigned j(0); j != theIntersectLipids[i].size(); ++j)
            {
              unsigned aMol(theIntersectLipids[i][j]+lipStartMol);
              if(theLattice[aMol*theAdjoinSize] == anID)
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
          thePopulatableMols.resize(0);
          for(unsigned i(0); i != theMolSize; ++i)
            {
              unsigned j(getMol(i)-vacStartMol);
              for(unsigned k(0); k != theIntersectLipids[j].size(); ++k)
                {
                  unsigned aMol(theIntersectLipids[j][k]+lipStartMol);
                  if(theLattice[aMol*theAdjoinSize] == theID)
                    {
                      thePopulatableMols.push_back(aMol);
                    }
                }
            }
          return thePopulatableMols.size();
        }
      return theMolSize;
    }
  unsigned getRandomPopulatableMol()
    {
      unsigned aMol;
      if(isMultiscale)
        {
          unsigned index(0);
          do
            {
              index = gsl_rng_uniform_int(theRng, thePopulatableMols.size());
            }
          while(theLattice[thePopulatableMols[index]*theAdjoinSize] != theID);
          aMol =  thePopulatableMols[index];
        }
      else
        {
          aMol = getRandomMol();
          while(theLattice[aMol*theAdjoinSize] != theID)
            {
              aMol = getRandomMol();
            }
        }
      return aMol;
    }
  unsigned getPopulatableMol(unsigned index)
    {
      if(isMultiscale)
        {
          return thePopulatableMols[index];
        }
      return getMol(index);
    }
  void setMultiscaleVacantSpecies(Species* aSpecies)
    {
      theMultiscaleVacantSpecies = aSpecies;
    }
  //Can aVoxel be populated by this species:
  bool isPopulatable(unsigned aMol)
    {
      if(isMultiscale)
        {
          if(isIntersectMultiscale(aMol))
            {
              return false;
            }
        }
      else if(theLattice[aMol*theAdjoinSize] != theVacantID)
        {
          return false;
        }
      return true;
    }
  //Can aVoxel of this species replaced by aSpecies:
  bool isReplaceable(unsigned aMol, Species* aSpecies)
    {
      if(getComp() != aSpecies->getComp() &&
         theID != aSpecies->getVacantID())
        {
          return false;
        }
      if(aSpecies->getIsMultiscale())
        {
          if(aSpecies->isIntersectMultiscale(aMol))
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
  const unsigned theID;
  unsigned lipStartMol;
  unsigned theAdjoinSize;
  unsigned theCollision;
  unsigned theDimension;
  unsigned theInitMolSize;
  unsigned theMolSize;
  unsigned theNullMol;
  unsigned theNullID;
  unsigned theSpeciesSize;
  unsigned vacStartMol;
  int thePolymerDirectionality;
  unsigned theVacantID;
  double D;
  double theDiffuseRadius;
  double theDiffusionInterval;
  double theMolRadius;
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
  std::vector<unsigned> theMols;
  std::vector<unsigned> theMultiscaleBindIDs;
  std::vector<unsigned> theMultiscaleIDs;
  std::vector<unsigned> theMultiscaleUnbindIDs;
  std::vector<unsigned> thePopulatableMols;
  std::vector<Tag> theTags;
  std::vector<double> theBendAngles;
  std::vector<double> theReactionProbabilities;
  std::vector<unsigned>* theCompMols;
  std::vector<Species*> theDiffusionInfluencedReactantPairs;
  std::vector<Species*> theTaggedSpeciesList;
  std::vector<Species*> theTagSpeciesList;
  std::vector<DiffusionInfluencedReactionProcess*> 
    theDiffusionInfluencedReactions;
  std::vector<SpatiocyteNextReactionProcess*> theInterruptedProcesses;
  std::vector<Origin> theMolOrigins;
  std::vector<unsigned> theRands;
  std::vector<unsigned> theTarMols;
  std::vector<VoxelInfo>& theInfo;
  std::vector<unsigned>& theLattice;
  std::vector<std::vector<unsigned> > theIntersectLipids;
};


#endif /* __SpatiocyteSpecies_hpp */

