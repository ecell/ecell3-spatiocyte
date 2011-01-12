//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of E-Cell Simulation Environment package
//
//                Copyright (C) 2006-2009 Keio University
//
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
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

#include <Variable.hpp>
#include "SpatiocyteCommon.hpp"
#include "SpatiocyteStepper.hpp"
#include "DiffusionInfluencedReactionProcess.hpp"
#include "PolymerFragmentationProcess.hpp"
#include "MoleculePopulateProcess.hpp"

// The size of Voxel must be 128 bytes to avoid cacheline splits
// The Core 2 has 64-byte cacheline
static double getDistance(Point* aSourcePoint, Point* aDestPoint)
{
  return sqrt(pow(aDestPoint->x-aSourcePoint->x, 2)+
              pow(aDestPoint->y-aSourcePoint->y, 2)+
              pow(aDestPoint->z-aSourcePoint->z, 2));
}

class Species
{
public:
  Species(SpatiocyteStepper* aStepper, Variable* aVariable, int anID, 
          int anInitMoleculeSize, const gsl_rng* aRng):
    isCentered(false),
    isDiffusing(false),
    isGaussianPopulation(false),
    isInContact(false),
    isLipid(false),
    isPolymer(false),
    isStatic(true),
    isSubunitInitialized(false),
    isVacant(false),
    isVolume(false),
    theID(anID),
    theInitMoleculeSize(anInitMoleculeSize),
    theMoleculeSize(0),
    D(0),
    theDiffusionInterval(libecs::INF),
    theWalkProbability(1),
    theRng(aRng),
    thePopulateProcess(NULL),
    theStepper(aStepper),
    theVariable(aVariable) {}
  ~Species() {}
  void initialize(int speciesSize)
    {
      theReactionProbabilities.resize(speciesSize);
      theDiffusionInfluencedReactions.resize(speciesSize);
      for(int i(0); i != speciesSize; ++ i)
        {
          theDiffusionInfluencedReactions[i] = NULL;
          theReactionProbabilities[i] = 0;
        }
    }
  void setDiffusionInfluencedReaction(DiffusionInfluencedReactionProcess*
                                      aReaction, int anID, double aProbability)
    {
      theDiffusionInfluencedReactions[anID] = aReaction;
      theReactionProbabilities[anID] = aProbability;
    }
  void setDiffusionInfluencedReactantPair(Species* aSpecies)
    {
      theDiffusionInfluencedReactantPairs.push_back(aSpecies);
    }
  void setPopulateProcess(MoleculePopulateProcess* aProcess, double aDist)
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
  void populateCompartmentGaussian()
    {
      if(thePopulateProcess)
        {
          thePopulateProcess->populateGaussian(this);
        }
      else if(theMoleculeSize)
        {
          std::cout << "Species:" << theVariable->getFullID().asString() <<
            " not MoleculePopulated." << std::endl;
        }
    }
  void populateCompartmentUniform(unsigned int voxelIDs[], unsigned int* aCount)
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
  void populateCompartmentUniformSparse()
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
  Variable* getVariable()
    {
      return theVariable;
    }
  std::vector<unsigned int> getSourceCoords()
    {
      std::vector<unsigned int> aCoords;
      for(unsigned int i(0); i != theMoleculeSize; ++i)
        {
          std::vector<Voxel*>& aSourceVoxels(theMolecules[i]->subunit->sourceVoxels);
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
          std::vector<Voxel*>& aTargetVoxels(theMolecules[i]->subunit->targetVoxels);
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
          std::vector<Voxel*>& aSharedLipids(theMolecules[i]->subunit->sharedLipids);
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
  const unsigned int size() const
    {
      return theMoleculeSize;
    }
  unsigned int getCoord(int anIndex)
    {
      return theMolecules[anIndex]->coord;
    }
  Point getPoint(int anIndex)
    {
      return theMolecules[anIndex]->subunit->subunitPoint;
    }
  Voxel* getMolecule(int anIndex)
    {
      return theMolecules[anIndex];
    }
  const unsigned short getID() const
    {
      return theID;
    }
  const double getMeanSquaredDisplacement()
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
                                                 isVolume,
                                                 &theMoleculeOrigins[i]));
          double aDistance(getDistance(&theMoleculeOrigins[i].point,
                                       &aCurrentPoint));
          aDisplacement += aDistance*aDistance;
        }
      return
        aDisplacement*pow(theStepper->getVoxelRadius()*2, 2)/theMoleculeSize;
    }
  void setIsLipid()
    {
      isLipid = true;
      isStatic = false;
    }
  void setIsSubunitInitialized()
    {
      isSubunitInitialized = true;
    }
  void setIsVacant()
    {
      isVacant = true;
      isStatic = false;
    }
  void setIsVolume()
    {
      isVolume = true;
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
      getVariable()->setValue(theMoleculeSize);
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
          isStatic = false;
        }
    }
  const double getDiffusionCoefficient() const
    {
      return D;
    }
  const double getWalkProbability() const
    {
      return theWalkProbability;
    }
  const bool getIsVolume() const
    {
      return isVolume;
    }
  const bool getIsStatic() const
    {
      return isStatic;
    }
  const bool getIsPolymer() const
    {
      return isPolymer;
    }
  const bool getIsSubunitInitialized() const
    {
      return isSubunitInitialized;
    }
  const bool getIsDiffusing() const
    {
      return isDiffusing;
    }
  const bool getIsLipid() const
    {
      return isLipid;
    }
  const bool getIsVacant() const
    {
      return isVacant;
    }
  const bool getIsInContact() const
    {
      return isInContact;
    }
  const bool getIsCentered() const
    {
      return isCentered;
    }
  const bool getIsPopulated() const
    {
      return theMoleculeSize == theInitMoleculeSize;
    }
  const double getDiffusionInterval() const
    {
      return theDiffusionInterval;
    }
  void surfaceWalkCollide()
    {
      int vacantID(theCompartment->vacantID);
      const int r(gsl_rng_uniform_int(theRng, 6));
      Voxel* target;
      for(unsigned int i(0); i < theMoleculeSize; ++i)
        {
          Voxel* source(theMolecules[i]);
          const unsigned short size(source->adjoiningSize);
          if(size > 6 || size <= r)
            { 
              if(size)
                {
                  target = source->adjoiningVoxels[
                    gsl_rng_uniform_int(theRng, size)];
                }
              else
                {
                  continue;
                }
            }
          else
            {
              target = source->adjoiningVoxels[r];
            }
          if(source == target)
            {
              std::cout << "error surface" << std::endl;
            }
          if(target->id == vacantID)
            {
              if(gsl_rng_uniform(theRng) < theWalkProbability)
                {
                  target->id = theID;
                  source->id = vacantID;
                  theMolecules[i] = target;
                }
            }
          else if(theDiffusionInfluencedReactions[target->id] != NULL)
            {
              //If it meets the reaction probability:
              if(gsl_rng_uniform(theRng) < theReactionProbabilities[target->id])
                { 
                  Species* targetSpecies(theStepper->id2species(target->id));
                  DiffusionInfluencedReactionProcess* aReaction(
                             theDiffusionInfluencedReactions[target->id]);
                  if(aReaction->react(source, &target))
                    {
                      //Soft remove the source molecule, i.e.,
                      //keep the id intact:
                      theMolecules[i--] = theMolecules[--theMoleculeSize];
                      theVariable->setValue(theMoleculeSize);
                      for(unsigned int j(0);
                          j != theInterruptedProcesses.size(); ++j)
                        {
                          theInterruptedProcesses[j
                            ]->removeSubstrateInterrupt(this, source);
                        }
                      //Soft remove the target molecule:
                      targetSpecies->softRemoveMolecule(target);
                      aReaction->finalizeReaction();
                    }
                }
            }
        }
    }
  void volumeWalkCollide()
    {
      int vacantID(theCompartment->vacantID);
      const int r(gsl_rng_uniform_int(theRng, ADJOINING_VOXEL_SIZE));
      for(unsigned int i(0); i < theMoleculeSize; ++i)
        {
          Voxel* source(theMolecules[i]);
          Voxel* target(source->adjoiningVoxels[r]);
          if(source == target)
            {
              //We cannot have a volume/surface molecule colliding with itself
              //because we need to ensure homodimerization reaction is
              //correct.
              std::cout << "error volume" << std::endl;
            }
          //walk:
          if(target->id == vacantID)
            {
              if(gsl_rng_uniform(theRng) < theWalkProbability)
                {
                  target->id = theID;
                  source->id = vacantID;
                  theMolecules[i] = target;
                }
            }
          //collide;
          else if(theDiffusionInfluencedReactions[target->id] != NULL)
            {
              //If it meets the reaction probability:
              if(gsl_rng_uniform(theRng) < theReactionProbabilities[target->id])
                { 
                  Species* targetSpecies(theStepper->id2species(target->id));
                  DiffusionInfluencedReactionProcess* aReaction(
                             theDiffusionInfluencedReactions[target->id]);
                  if(aReaction->react(source, &target))
                    {
                      //Soft remove the source molecule, i.e.,
                      //keep the id intact:
                      theMolecules[i--] = theMolecules[--theMoleculeSize];
                      theVariable->setValue(theMoleculeSize);
                      for(unsigned int j(0);
                          j != theInterruptedProcesses.size(); ++j)
                        {
                          theInterruptedProcesses[j
                            ]->removeSubstrateInterrupt(this, source);
                        }
                      //Soft remove the target molecule:
                      targetSpecies->softRemoveMolecule(target);
                      aReaction->finalizeReaction();
                    }
                }
            }
        }
    }
  void surfaceWalk()
    {
      int vacantID(theCompartment->vacantID);
      const int r(gsl_rng_uniform_int(theRng, 6));
      Voxel* target;
      for(unsigned int i(0); i < theMoleculeSize; ++i)
        {
          Voxel* source(theMolecules[i]);
          const unsigned short size(source->adjoiningSize);
          if(size > 6 || size <= r)
            { 
              if(size)
                {
                  target = source->adjoiningVoxels[
                    gsl_rng_uniform_int(theRng, size)];
                }
              else
                {
                  continue;
                }
            }
          else
            {
              target = source->adjoiningVoxels[r];
            }
          if(source == target)
            {
              std::cout << "error surface" << std::endl;
            }
          if(target->id == vacantID)
            {
              target->id = theID;
              source->id = vacantID;
              theMolecules[i] = target;
            }
          else if(theDiffusionInfluencedReactions[target->id] != NULL)
            {
              //If it meets the reaction probability:
              if(gsl_rng_uniform(theRng) < theReactionProbabilities[target->id])
                { 
                  Species* targetSpecies(theStepper->id2species(target->id));
                  DiffusionInfluencedReactionProcess* aReaction(
                             theDiffusionInfluencedReactions[target->id]);
                  if(aReaction->react(source, &target))
                    {
                      //Soft remove the source molecule, i.e.,
                      //keep the id intact:
                      theMolecules[i--] = theMolecules[--theMoleculeSize];
                      theVariable->setValue(theMoleculeSize);
                      for(unsigned int j(0);
                          j != theInterruptedProcesses.size(); ++j)
                        {
                          theInterruptedProcesses[j
                            ]->removeSubstrateInterrupt(this, source);
                        }
                      //Soft remove the target molecule:
                      targetSpecies->softRemoveMolecule(target);
                      aReaction->finalizeReaction();
                    }
                }
            }
        }
    }
  void volumeWalk()
    {
      int vacantID(theCompartment->vacantID);
      const int r(gsl_rng_uniform_int(theRng, ADJOINING_VOXEL_SIZE));
      for(unsigned int i(0); i < theMoleculeSize; ++i)
        {
          Voxel* source(theMolecules[i]);
          Voxel* target(source->adjoiningVoxels[r]);
          if(source == target)
            {
              //We cannot have a volume/surface molecule colliding with itself
              //because we need to ensure homodimerization reaction is
              //correct.
              std::cout << "error volume" << std::endl;
            }
          if(target->id == vacantID)
            {
              target->id = theID;
              source->id = vacantID;
              theMolecules[i] = target;
            }
          else if(theDiffusionInfluencedReactions[target->id] != NULL)
            {
              //If it meets the reaction probability:
              if(gsl_rng_uniform(theRng) < theReactionProbabilities[target->id])
                { 
                  Species* targetSpecies(theStepper->id2species(target->id));
                  DiffusionInfluencedReactionProcess* aReaction(
                             theDiffusionInfluencedReactions[target->id]);
                  if(aReaction->react(source, &target))
                    {
                      //Soft remove the source molecule, i.e.,
                      //keep the id intact:
                      theMolecules[i--] = theMolecules[--theMoleculeSize];
                      theVariable->setValue(theMoleculeSize);
                      for(unsigned int j(0);
                          j != theInterruptedProcesses.size(); ++j)
                        {
                          theInterruptedProcesses[j
                            ]->removeSubstrateInterrupt(this, source);
                        }
                      //Soft remove the target molecule:
                      targetSpecies->softRemoveMolecule(target);
                      aReaction->finalizeReaction();
                    }
                }
            }
        }
    }
  void setCompartment(Compartment* aCompartment)
    {
      theCompartment = aCompartment;
    }
  Compartment* getCompartment()
    {
      return theCompartment;
    }
  void setVariable(Variable* aVariable)
    {
      theVariable = aVariable;
    }
  void addSimpleMolecule(Voxel* aMolecule)
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
  void addMolecule(Voxel* aMolecule)
    {
      if(!isVacant && !isLipid)
        {
          ++theMoleculeSize;
          aMolecule->id = theID;
          if(theMoleculeSize > theMolecules.size())
            {
              theMolecules.push_back(aMolecule);
            }
          else
            {
              theMolecules[theMoleculeSize-1] = aMolecule;
            }
          theVariable->setValue(theMoleculeSize);
          for(unsigned int i(0); i != theInterruptedProcesses.size(); ++i)
            {
              theInterruptedProcesses[i]->addSubstrateInterrupt(this,
                                                                aMolecule);
            }
        }
    }
  void softRemoveMolecule(Voxel* aMolecule)
    {
      for(unsigned int i(0); i < theMoleculeSize; ++i)
        {
          if(theMolecules[i] == aMolecule)
            {
              theMolecules[i] = theMolecules[--theMoleculeSize];
              theVariable->setValue(theMoleculeSize);
              for(unsigned int i(0); i != theInterruptedProcesses.size(); ++i)
                {
                  theInterruptedProcesses[i]->removeSubstrateInterrupt(this,
                                                                     aMolecule);
                }
              return;
            }
        }
      if(!isLipid)
        {
          if(aMolecule != aMolecule->subunit->voxel)
            {
              std::cout << "it is a shared molecule" << std::endl;
            }
          std::cout << "error in soft removing molecule, couldn't find the"
           " specified molecule to be removed, molecule size:" << 
            theMoleculeSize << " value:" << theVariable->getValue() <<
            getVariable()->getFullID().asString() << std::endl;
          std::cout << "species:" << theStepper->id2species(aMolecule->id)->getVariable()->getFullID().asString() << std::endl;
        }
    }
  void removeMolecule(Voxel* aMolecule)
    {
      for(unsigned int i(0); i < theMoleculeSize; ++i)
        {
          if(theMolecules[i] == aMolecule)
            {
              aMolecule->id = theCompartment->vacantID;
              theMolecules[i] = theMolecules[--theMoleculeSize];
              theVariable->setValue(theMoleculeSize);
              for(unsigned int i(0); i != theInterruptedProcesses.size(); ++i)
                {
                  theInterruptedProcesses[i]->removeSubstrateInterrupt(this,
                                                                     aMolecule);
                }
              return;
            }
        }
      if(!isLipid)
        {
          std::cout << "error in removing molecule, couldn't find the specified" <<
            " molecule to be removed, molecule size:" << 
            theMoleculeSize << " value:" << theVariable->getValue() <<
            getVariable()->getFullID().asString() << std::endl;
        }
    }
  void removeMolecules()
    {
      for(unsigned int i(0); i < theMoleculeSize; ++i)
        {
          theMolecules[i]->id = theCompartment->vacantID;
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
          if(theStepper->isBoundaryCoord(theMolecules[i]->coord, isVolume))
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
                                                isVolume,
                                                &anOrigin));
          if(periodicVoxel != NULL && 
             periodicVoxel->id == theCompartment->vacantID)
            {
              theMolecules[i]->id = theCompartment->vacantID;
              theMolecules[i] = periodicVoxel;
              theMolecules[i]->id = theID;
              theMoleculeOrigins[i] = anOrigin;
            }
        }
    }
  const int getVacantID() const
    {
      return theCompartment->vacantID;
    }
  const std::vector<double>& getBendAngles() const
    {
      return theBendAngles;
    }
  const Point& getWestPoint() const
    {
      return theCompartment->westPoint;
    }
  const Point& getEastPoint() const
    {
      return theCompartment->eastPoint;
    }
  const double getRadius() const
    {
      return theCompartment->lengthY/2;
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
          std::cout << "size:" << theVariable->getValue() << std::endl;
        }
      return theMolecules[gsl_rng_uniform_int(theRng, theMoleculeSize)];
    }
  void addInterruptedProcess(SpatiocyteProcess* aProcess)
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
      std::vector<Voxel*> compartmentVoxels;
      if(theStepper->getSearchVacant())
        { 
          for(unsigned int i(0); i != ADJOINING_VOXEL_SIZE; ++i)
            {
              Voxel* aVoxel(source->adjoiningVoxels[i]);
              if(aVoxel->id == theCompartment->vacantID)
                {
                  compartmentVoxels.push_back(aVoxel);
                }
            }
        }
      else
        {
          for(unsigned int i(0); i != ADJOINING_VOXEL_SIZE; ++i)
            {
              Voxel* aVoxel(source->adjoiningVoxels[i]);
              if(theStepper->id2compartment(aVoxel->id) == theCompartment)
                {
                  compartmentVoxels.push_back(aVoxel);
                }
            }
        }
      return getRandomVacantVoxel(&compartmentVoxels);
    } 
  Voxel* getRandomAdjoiningVoxel(Voxel* source, Voxel* target)
    {
      std::vector<Voxel*> compartmentVoxels;
      if(theStepper->getSearchVacant())
        { 
          for(unsigned int i(0); i != ADJOINING_VOXEL_SIZE; ++i)
            {
              Voxel* aVoxel(source->adjoiningVoxels[i]);
              if(aVoxel->id == theCompartment->vacantID &&
                 aVoxel != target)
                {
                  compartmentVoxels.push_back(aVoxel);
                }
            }
        }
      else
        {
          for(unsigned int i(0); i != ADJOINING_VOXEL_SIZE; ++i)
            {
              Voxel* aVoxel(source->adjoiningVoxels[i]);
              if(theStepper->id2compartment(aVoxel->id) == theCompartment &&
                 aVoxel != target)
                {
                  compartmentVoxels.push_back(aVoxel);
                }
            }
        }
      return getRandomVacantVoxel(&compartmentVoxels);
    }
  Voxel* getRandomAdjoiningVoxel(Voxel* source, Voxel* targetA, Voxel* targetB)
    {
      std::vector<Voxel*> compartmentVoxels;
      if(theStepper->getSearchVacant())
        { 
          for(unsigned int i(0); i != ADJOINING_VOXEL_SIZE; ++i)
            {
              Voxel* aVoxel(source->adjoiningVoxels[i]);
              if(aVoxel->id == theCompartment->vacantID &&
                 aVoxel != targetA && aVoxel != targetB)
                {
                  compartmentVoxels.push_back(aVoxel);
                }
            }
        }
      else
        {
          for(unsigned int i(0); i != ADJOINING_VOXEL_SIZE; ++i)
            {
              Voxel* aVoxel(source->adjoiningVoxels[i]);
              if(theStepper->id2compartment(aVoxel->id) == theCompartment &&
                 aVoxel != targetA && aVoxel != targetB)
                {
                  compartmentVoxels.push_back(aVoxel);
                }
            }
        }
      return getRandomVacantVoxel(&compartmentVoxels);
    }
  Voxel* getRandomVacantVoxel(std::vector<Voxel*>* voxels)
    {
      if(voxels->size())
        {
          const int r(gsl_rng_uniform_int(theRng, voxels->size())); 
          Voxel* aVoxel((*voxels)[r]);
          if(aVoxel->id == theCompartment->vacantID)
            {
              return aVoxel;
            }
        }
      return NULL;
    }
  Voxel* getRandomCompartmentVoxel()
    {
      int aSize(theCompartment->coords.size());
      int r(gsl_rng_uniform_int(theRng, aSize));
      if(theStepper->getSearchVacant())
        {
          for(int i(r); i != aSize; ++i)
            {
              Voxel* aVoxel(theStepper->coord2voxel(theCompartment->coords[i]));
              if(aVoxel->id == theCompartment->vacantID)
                {
                  return aVoxel;
                }
            }
          for(int i(0); i != r; ++i)
            {
              Voxel* aVoxel(theStepper->coord2voxel(theCompartment->coords[i]));
              if(aVoxel->id == theCompartment->vacantID)
                {
                  return aVoxel;
                }
            }
        }
      else
        {
          Voxel* aVoxel(theStepper->coord2voxel(theCompartment->coords[r]));
          if(aVoxel->id == theCompartment->vacantID)
            {
              return aVoxel;
            }
        }
      return NULL;
    }
private:
  bool isCentered;
  bool isDiffusing;
  bool isGaussianPopulation;
  bool isInContact;
  bool isLipid;
  bool isPolymer;
  bool isStatic;
  bool isSubunitInitialized;
  bool isVacant;
  bool isVolume;
  const unsigned short theID;
  const unsigned int theInitMoleculeSize;
  unsigned int theMoleculeSize;
  int thePolymerDirectionality;
  double D;
  double theDiffusionInterval;
  double theWalkProbability;
  const gsl_rng* theRng;
  Compartment* theCompartment;
  MoleculePopulateProcess* thePopulateProcess;
  SpatiocyteStepper* theStepper;
  Variable* theVariable;
  std::vector<double> theBendAngles;
  std::vector<double> theReactionProbabilities;
  std::vector<Voxel*> theMolecules;
  std::vector<Species*> theDiffusionInfluencedReactantPairs;
  std::vector<DiffusionInfluencedReactionProcess*> theDiffusionInfluencedReactions;
  std::vector<SpatiocyteProcess*> theInterruptedProcesses;
  std::vector<Origin> theMoleculeOrigins;
};






#endif /* __SpatiocyteSpecies_hpp */
