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


#ifndef __SpatiocyteProcess_hpp
#define __SpatiocyteProcess_hpp

#include <Process.hpp>
#include <Model.hpp>
#include "SpatiocyteCommon.hpp"
#include "SpatiocyteStepper.hpp"
#include "SpatiocyteProcessInterface.hpp"

class ReactionProcess;

LIBECS_DM_CLASS_EXTRA_1(SpatiocyteProcess, Process, virtual SpatiocyteProcessInterface)
{ 
public:
  LIBECS_DM_OBJECT(SpatiocyteProcess, Process)
    {
      INHERIT_PROPERTIES(Process);
    }
  SpatiocyteProcess():
    isExternInterrupted(false),
    isInitialized(false),
    isPriorityQueued(false),
    theInterval(libecs::INF),
    theTime(libecs::INF) {}
  virtual ~SpatiocyteProcess() {}
  virtual void fire() {}
  virtual void initializeFirst() {}
  virtual void initializeSecond()
    {
      theSpecies = theSpatiocyteStepper->getSpecies();
    }
  virtual void initializeThird() {}
  virtual void initializeFourth() {}
  virtual void initializeFifth() {}
  virtual void initializeLastOnce() {}
  virtual void printParameters() {}
  virtual void updateResizedLattice() {}
  virtual void prepreinitialize()
    {
      //Process::prepreinitialize();
    }
  virtual void initialize()
    {
      if(isInitialized)
        {
          return;
        }
      isInitialized = true;
      Process::initialize();
      theSpatiocyteStepper = dynamic_cast<SpatiocyteStepper*>(getStepper());
      theSortedVariableReferences.resize(theVariableReferenceVector.size());
      for(VariableReferenceVector::iterator
          i(theVariableReferenceVector.begin());
          i != theVariableReferenceVector.end(); ++i)
        {
          theSortedVariableReferences[(*i).getEllipsisNumber()] = *i; 
        }
      for(VariableReferenceVector::iterator
          i(theSortedVariableReferences.begin());
          i != theSortedVariableReferences.end(); ++i)
        {
          Variable* aVariable((*i).getVariable()); 
          // if Name is HD, it is a homogeneously distributed species.
          // only create Species for immobile or diffusing species
          if(aVariable->getName() != "HD")
            {
              Species* aSpecies(theSpatiocyteStepper->addSpecies(aVariable));
              theProcessSpecies.push_back(aSpecies);
            }
          if((*i).getCoefficient() > 0)
            {
              thePositiveVariableReferences.push_back(*i);
            }
          else if((*i).getCoefficient() < 0)
            {
              theNegativeVariableReferences.push_back(*i);
            }
          else
            {
              theZeroVariableReferences.push_back(*i);
            }
        }
    }
  void requeue()
    {
      theTime += getInterval(); // do this only for the Processes in Q
      thePriorityQueue->moveTop(); // do this only for the Processes in Q
    }
  virtual bool isInterrupted(ReactionProcess* aProcess)
    {
      return false;
    }
  virtual double getInterval()
    {
      return theInterval;
    }
  virtual void setPriorityQueue(ProcessPriorityQueue* aPriorityQueue)
    {
      thePriorityQueue = aPriorityQueue;
    }
  virtual void setLatticeProperties(std::vector<Voxel>* aLattice,
                                    unsigned anAdjoiningCoordSize,
                                    unsigned aNullCoord, unsigned aNullID,
                                    RandomLib::Random* aRng)
    {
      theRng = aRng;
      theLattice = aLattice;
      theAdjoiningCoordSize = anAdjoiningCoordSize;
      theNullCoord = aNullCoord;
      theNullID = aNullID;
    }
  double getTime() const
    {
      return theTime;
    }
  virtual int getPriority() const
    {
      return Process::getPriority();
    }
  virtual void setTime(double aTime)
    {
      theTime = aTime;
    }
  virtual void setQueueID(ProcessID anID)
    {
      theQueueID = anID;
    }
  Species* id2species(unsigned short id)
    {
      return theSpatiocyteStepper->id2species(id);
    }
  virtual unsigned getLatticeResizeCoord(unsigned)
    {
      return 0;
    }
  virtual void addSubstrateInterrupt(Species* aSpecies, Voxel* aMolecule) {}
  virtual void removeSubstrateInterrupt(Species* aSpecies, Voxel* aMolecule) {}
  virtual void substrateValueChanged(double aCurrentTime)
    {
      const double anOldTime(theTime);
      theTime = aCurrentTime + getInterval();
      if(theTime >= anOldTime)
        {
          thePriorityQueue->moveDown(theQueueID);
        }
      else if(theTime < anOldTime)
        {
          thePriorityQueue->moveUp(theQueueID);
        }          
    }
  virtual bool getIsExternInterrupted()
    {
      return isExternInterrupted;
    }
  virtual bool getIsPriorityQueued()
    {
      return isPriorityQueued;
    }
  Variable* createVariable(String anID)
    {
      String anEntityType("Variable");
      SystemPath aSystemPath(getSuperSystem()->getSystemPath());
      aSystemPath.push_back(getSuperSystem()->getID());
      FullID aFullID(anEntityType, aSystemPath, anID);
      Variable* aVariable(reinterpret_cast<Variable*>(
                              getModel()->createEntity("Variable", aFullID)));
      aVariable->setValue(0);
      return aVariable;
    }
  String getIDString()
    {
      return getFullID().asString();
    }
protected:
  String getIDString(Voxel*) const;
  String getIDString(Species*) const;
  String getIDString(Variable*) const;
  String getIDString(Comp*) const;
  String getIDString(unsigned) const;
protected:
  bool isExternInterrupted;
  bool isInitialized;
  bool isPriorityQueued;
  unsigned theAdjoiningCoordSize;
  unsigned theNullCoord;
  unsigned theNullID;
  double theInterval;
  double theTime;
  ProcessID theQueueID;
  ProcessPriorityQueue* thePriorityQueue; 
  SpatiocyteStepper* theSpatiocyteStepper;
  std::vector<Species*> theSpecies;
  std::vector<Species*> theProcessSpecies;
  std::vector<Voxel>* theLattice;
  VariableReferenceVector thePositiveVariableReferences;
  VariableReferenceVector theNegativeVariableReferences;
  VariableReferenceVector theZeroVariableReferences;
  VariableReferenceVector theSortedVariableReferences;
  RandomLib::Random* theRng;
};

#endif /* __SpatiocyteProcess_hpp */
