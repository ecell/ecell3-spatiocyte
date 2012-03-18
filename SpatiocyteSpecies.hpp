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

class Species
{
public:
  Species(SpatiocyteStepper* aStepper, Variable* aVariable, int anID, 
          int anInitMoleculeSize, const gsl_rng* aRng):
      isVacant(false),
      isVolume(false),
    isCentered(false),
    isDiffusing(false),
    isGaussianPopulation(false),
    isInContact(false),
    isPolymer(false),
    isStatic(true),
    isSubunitInitialized(false),
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
  void initialize(int speciesSize, int anAdjoiningVoxelSize)
    {
      theAdjoiningVoxelSize = anAdjoiningVoxelSize;
      theReactionProbabilities.resize(speciesSize);
      theDiffusionInfluencedReactions.resize(speciesSize);
      for(int i(0); i != speciesSize; ++ i)
        {
          theDiffusionInfluencedReactions[i] = NULL;
          theReactionProbabilities[i] = 0;
        }
      if(theComp)
        {
          setVacantSpecies(theStepper->id2species(theComp->vacantID));
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
  void populateCompGaussian()
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
  Variable* getVariable() const
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
      return theMolecules[anIndex]->subunit->subunitPoint;
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
                                                 getIsVolume(),
                                                 &theMoleculeOrigins[i]));
          double aDistance(getDistance(&theMoleculeOrigins[i].point,
                                       &aCurrentPoint));
          aDisplacement += aDistance*aDistance;
        }
      return
        aDisplacement*pow(theStepper->getVoxelRadius()*2, 2)/theMoleculeSize;
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
  double getDiffusionCoefficient() const
    {
      return D;
    }
  double getWalkProbability() const
    {
      return theWalkProbability;
    }
  bool getIsVolume() const
    {
        return isVolume;
    }
  bool getIsStatic() const
    {
      return isStatic;
    }
  bool getIsPolymer() const
    {
      return isPolymer;
    }
  bool getIsSubunitInitialized() const
    {
      return isSubunitInitialized;
    }
  bool getIsDiffusing() const
    {
      return isDiffusing;
    }
  bool getIsVacant() const
    {
      return isVacant;
    }
  bool getIsLipid() const
    {
      return (isVacant && theComp->dimension == 2);
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
  void surfaceWalkCollide()
    {
      Voxel* target;
      for(unsigned int i(0); i < theMoleculeSize; ++i)
        {
          Voxel* source(theMolecules[i]);
          const unsigned short size(source->adjoiningSize);
          target = source->adjoiningVoxels[gsl_rng_uniform_int(theRng, size)];
          if(source == target)
            {
              std::cout << "error surface" << std::endl;
            }
          if(target->id == theVacantID)
            {
              if(gsl_rng_uniform(theRng) < theWalkProbability)
                {
                  theVacantSpecies->softReplaceMolecule(target, source);
                  target->id = theID;
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
      for(unsigned int i(0); i < theMoleculeSize; ++i)
        {
          Voxel* source(theMolecules[i]); 
          const int r(gsl_rng_uniform_int(theRng, theAdjoiningVoxelSize));
          Voxel* target(source->adjoiningVoxels[r]);
          if(source == target)
            {
              //We cannot have a volume/surface molecule colliding with itself
              //because we need to ensure homodimerization reaction is
              //correct.
              std::cout << "error volume" << std::endl;
            }
          //walk:
          if(target->id == theVacantID)
            {
              if(gsl_rng_uniform(theRng) < theWalkProbability)
                {
                  theVacantSpecies->softReplaceMolecule(target, source);
                  target->id = theID;
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
                  DiffusionInfluencedReactionProcessInterface* aReaction(
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
      Voxel* target;
      for(unsigned int i(0); i < theMoleculeSize; ++i)
        {
          Voxel* source(theMolecules[i]);
          const unsigned short size(source->adjoiningSize); 
          target = source->adjoiningVoxels[gsl_rng_uniform_int(theRng, size)];
          if(source == target)
            {
              std::cout << "error surface" << std::endl;
            }
          if(target->id == theVacantID)
            {
              theVacantSpecies->softReplaceMolecule(target, source);
              target->id = theID;
              theMolecules[i] = target;
            }
          else if(theDiffusionInfluencedReactions[target->id] != NULL)
            {
              //If it meets the reaction probability:
              if(gsl_rng_uniform(theRng) < theReactionProbabilities[target->id])
                { 
                  Species* targetSpecies(theStepper->id2species(target->id));
                  DiffusionInfluencedReactionProcessInterface* aReaction(
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
      for(unsigned int i(0); i < theMoleculeSize; ++i)
        {
          Voxel* source(theMolecules[i]);
          const int r(gsl_rng_uniform_int(theRng, theAdjoiningVoxelSize));
          Voxel* target(source->adjoiningVoxels[r]);
          if(source == target)
            {
              //We cannot have a volume/surface molecule colliding with itself
              //because we need to ensure homodimerization reaction is
              //correct.
              std::cout << "error volume" << std::endl;
            }
          if(target->id == theVacantID)
            {
              theVacantSpecies->softReplaceMolecule(target, source);
              target->id = theID;
              theMolecules[i] = target;
            }
          else if(theDiffusionInfluencedReactions[target->id] != NULL)
            {
              //If it meets the reaction probability:
              if(gsl_rng_uniform(theRng) < theReactionProbabilities[target->id])
                { 
                  Species* targetSpecies(theStepper->id2species(target->id));
                  DiffusionInfluencedReactionProcessInterface* aReaction(
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
                          theInterruptedProcesses[j]->removeSubstrateInterrupt(this, source);
                        }
                      //Soft remove the target molecule:
                      targetSpecies->softRemoveMolecule(target);
                      aReaction->finalizeReaction();
                    }
                }
            }
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
      Species* aSpecies(theStepper->id2species(aMolecule->id));
      aMolecule->id = theID;
      if(!getIsVacant())
        {
          aSpecies->softRemoveMolecule(aMolecule);
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
          for(unsigned int i(0); i != theInterruptedProcesses.size(); ++i)
            {
              theInterruptedProcesses[i]->addSubstrateInterrupt(this,
                                                                aMolecule);
            }
        }
    }
  //it is soft replace because the id of the source molecule is not changed:
  void softReplaceMolecule(Voxel* aSource, Voxel* aTarget)
    {
      aTarget->id = theID;
      if(!getIsVacant())
        {
          for(unsigned int i(0); i < theMoleculeSize; ++i)
            {
              if(theMolecules[i] == aSource)
                {
                  theMolecules[i] = aTarget;
                  return;
                }
            }
        }
    }
  //it is soft remove because the id of the molecule is not changed:
  void softRemoveMolecule(Voxel* aMolecule)
    {
      if(!getIsVacant())
        {
          for(unsigned int i(0); i < theMoleculeSize; ++i)
            {
              if(theMolecules[i] == aMolecule)
                {
                  theMolecules[i] = theMolecules[--theMoleculeSize];
                  theVariable->setValue(theMoleculeSize);
                  for(unsigned int i(0); i != theInterruptedProcesses.size(); ++i)
                    {
                      theInterruptedProcesses[i]->removeSubstrateInterrupt(
                                                             this, aMolecule);
                    }
                  return;
                }
            }
        }
    }
  void removeMolecule(Voxel* aMolecule)
    {
      if(!getIsVacant())
        {
          for(unsigned int i(0); i < theMoleculeSize; ++i)
            {
              if(theMolecules[i] == aMolecule)
                {
                  theVacantSpecies->addMolecule(aMolecule);
                  theMolecules[i] = theMolecules[--theMoleculeSize];
                  theVariable->setValue(theMoleculeSize);
                  for(unsigned int i(0); 
                      i != theInterruptedProcesses.size(); ++i)
                    {
                      theInterruptedProcesses[i]->removeSubstrateInterrupt(
                                                             this, aMolecule);
                    }
                  return;
                }
            }
        }
    }
  //Used by the SpatiocyteStepper when resetting an interation, so must
  //clear the whole compartment:
  void removeMolecules()
    {
      for(unsigned int i(0); i < theMoleculeSize; ++i)
        {
          theMolecules[i]->id = theComp->vacantID;
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
                   theMolecules[i]->coord, getIsVolume()))
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
                                                getIsVolume(),
                                                &anOrigin));
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
  double getRadius() const
    {
      return theComp->lengthY/2;
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
          for(unsigned int i(0); i != theAdjoiningVoxelSize; ++i)
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
          for(unsigned int i(0); i != theAdjoiningVoxelSize; ++i)
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
  Voxel* getRandomAdjoiningVoxel(Voxel* source, Species* aVacantSpecies)
    {
      std::vector<Voxel*> CompVoxels;
      if(theStepper->getSearchVacant())
        { 
          for(unsigned int i(0); i != theAdjoiningVoxelSize; ++i)
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
          for(unsigned int i(0); i != theAdjoiningVoxelSize; ++i)
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
          for(unsigned int i(0); i != theAdjoiningVoxelSize; ++i)
            {
              Voxel* aVoxel(source->adjoiningVoxels[i]);
              if(aVoxel->id == theVacantID &&
                 aVoxel != target)
                {
                  CompVoxels.push_back(aVoxel);
                }
            }
        }
      else
        {
          for(unsigned int i(0); i != theAdjoiningVoxelSize; ++i)
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
          for(unsigned int i(0); i != theAdjoiningVoxelSize; ++i)
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
          for(unsigned int i(0); i != theAdjoiningVoxelSize; ++i)
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
      int aSize(theComp->coords.size());
      int r(gsl_rng_uniform_int(theRng, aSize));
      if(theStepper->getSearchVacant())
        {
          for(int i(r); i != aSize; ++i)
            {
              Voxel* aVoxel(theStepper->coord2voxel(theComp->coords[i]));
              if(aVoxel->id == theVacantID)
                {
                  return aVoxel;
                }
            }
          for(int i(0); i != r; ++i)
            {
              Voxel* aVoxel(theStepper->coord2voxel(theComp->coords[i]));
              if(aVoxel->id == theVacantID)
                {
                  return aVoxel;
                }
            }
        }
      else
        {
          Voxel* aVoxel(theStepper->coord2voxel(theComp->coords[r]));
          if(aVoxel->id == theVacantID)
            {
              return aVoxel;
            }
        }
      return NULL;
    }
  Voxel* getRandomAdjoiningCompVoxel(Comp* aComp)
    {
      int aSize(aComp->coords.size());
      int r(gsl_rng_uniform_int(theRng, aSize)); 
      Voxel* aVoxel(theStepper->coord2voxel(aComp->coords[r]));
      return getRandomAdjoiningVoxel(aVoxel);
    }
private:
  bool isVacant;
  bool isVolume;
  bool isCentered;
  bool isDiffusing;
  bool isGaussianPopulation;
  bool isInContact;
  bool isPolymer;
  bool isStatic;
  bool isSubunitInitialized;
  const unsigned short theID;
  const unsigned int theInitMoleculeSize;
  unsigned int theMoleculeSize;
  unsigned int theAdjoiningVoxelSize;
  int thePolymerDirectionality;
  int theVacantID;
  double D;
  double theDiffusionInterval;
  double theWalkProbability;
  const gsl_rng* theRng;
  Species* theVacantSpecies;
  Comp* theComp;
  MoleculePopulateProcessInterface* thePopulateProcess;
  SpatiocyteStepper* theStepper;
  Variable* theVariable;
  std::vector<double> theBendAngles;
  std::vector<double> theReactionProbabilities;
  std::vector<Voxel*> theMolecules;
  std::vector<Species*> theDiffusionInfluencedReactantPairs;
  std::vector<DiffusionInfluencedReactionProcessInterface*> theDiffusionInfluencedReactions;
  std::vector<SpatiocyteProcessInterface*> theInterruptedProcesses;
  std::vector<Origin> theMoleculeOrigins;
};






#endif /* __SpatiocyteSpecies_hpp */
