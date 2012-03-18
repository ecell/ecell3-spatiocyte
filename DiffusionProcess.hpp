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


#ifndef __DiffusionProcess_hpp
#define __DiffusionProcess_hpp

#include <sstream>
#include <MethodProxy.hpp>
#include "SpatiocyteProcess.hpp"
#include "SpatiocyteSpecies.hpp"

LIBECS_DM_CLASS(DiffusionProcess, SpatiocyteProcess)
{ 
  typedef void (DiffusionProcess::*WalkMethod)(void) const;
public:
  LIBECS_DM_OBJECT(DiffusionProcess, Process)
    {
      INHERIT_PROPERTIES(Process);
      PROPERTYSLOT_SET_GET(Real, D);
      PROPERTYSLOT_SET_GET(Real, P);
      PROPERTYSLOT_SET_GET(Real, WalkProbability);
    }
  DiffusionProcess():
    D(0),
    P(1),
    WalkProbability(1),
    theVacantSpecies(NULL),
    theWalkMethod(&DiffusionProcess::volumeWalk) {}
  virtual ~DiffusionProcess() {}
  SIMPLE_SET_GET_METHOD(Real, D);
  SIMPLE_SET_GET_METHOD(Real, P);
  SIMPLE_SET_GET_METHOD(Real, WalkProbability);
  virtual void initialize()
    {
      if(isInitialized)
        {
          return;
        }
      SpatiocyteProcess::initialize();
      for(VariableReferenceVector::iterator
          i(theVariableReferenceVector.begin());
          i != theVariableReferenceVector.end(); ++i)
        {
          Species* aSpecies(theSpatiocyteStepper->variable2species(
                                   (*i).getVariable())); 
          if(!(*i).getCoefficient())
            {
              theDiffusionSpecies.push_back(aSpecies);
              aSpecies->setDiffusionCoefficient(D);
            }
          else
            {
              theVacantSpecies = aSpecies;
            }
        }
      if(!theDiffusionSpecies.size())
        {
          THROW_EXCEPTION(ValueError, String(
                          getPropertyInterface().getClassName()) +
                          "[" + getFullID().asString() + 
                          "]: A DiffusionProcess requires at least one " +
                          "variable reference with zero coefficient."); 
        }
    }
  virtual void initializeSecond()
    {
      for(std::vector<Species*>::const_iterator
          i(theDiffusionSpecies.begin());
          i != theDiffusionSpecies.end(); ++i)
        {
          if(theVacantSpecies)
            {
              (*i)->setVacantSpecies(theVacantSpecies);
            }
        }
    }
  virtual void initializeThird()
    {
      Species* aSpecies(theDiffusionSpecies[0]);
      isVolume = aSpecies->getIsVolume();
      double rho(aSpecies->getMaxReactionProbability());
      if(D > 0)
        {
          for(std::vector<Species*>::const_iterator
              i(theDiffusionSpecies.begin());
              i != theDiffusionSpecies.end(); ++i)
            {
              if((*i)->getIsVolume() != isVolume)
                {
                  THROW_EXCEPTION(ValueError, String(
                                   getPropertyInterface().getClassName()) +
                                  "[" + getFullID().asString() + 
                                  "]: A DiffusionProcess can only execute" +
                                  " multiple species when they are all either" +
                                  " in a volume compartment or a surface" +
                                  " compartment, not both concurrently. " +
                                  getIDString(theDiffusionSpecies[0]) + " and " +
                                  getIDString(*i) + " belong to different" +
                                  " types of compartment.");
                }
              if(rho < (*i)->getMaxReactionProbability())
                {
                  if(rho > P)
                    {
                      THROW_EXCEPTION(ValueError, String(
                                       getPropertyInterface().getClassName()) + 
                                      "[" + getFullID().asString() + 
                                      "]: Create separate" +
                                      " DiffusionProcesses for " +
                                      getIDString(aSpecies) + " and " +
                                      getIDString(*i) + " since their" +
                                      " reaction probabilities are not the" +
                                      " same and the latter's reaction" +
                                      " probability is higher than P.");
                    }
                  aSpecies = *i;
                  rho = (*i)->getMaxReactionProbability();
                }
            }
        }
      if(rho > P)
        {
          WalkProbability = P/rho;
        }
      for(std::vector<Species*>::const_iterator i(theDiffusionSpecies.begin());
          i != theDiffusionSpecies.end(); ++i)
        {
          (*i)->rescaleReactionProbabilities(WalkProbability);
        }
      if(D > 0)
        {
          double r_v(theSpatiocyteStepper->getVoxelRadius());
          double lambda(2.0/3);
          if(!isVolume)
            {
              lambda = pow((2*sqrt(2)+4*sqrt(3)+3*sqrt(6)+sqrt(22))/
                           (6*sqrt(2)+4*sqrt(3)+3*sqrt(6)), 2);
            }
          theStepInterval = lambda*r_v*r_v*WalkProbability/D;
        }
      for(std::vector<Species*>::const_iterator i(theDiffusionSpecies.begin());
          i != theDiffusionSpecies.end(); ++i)
        {
          (*i)->setDiffusionInterval(theStepInterval);
        }
      if(isVolume)
        {
          if(WalkProbability == 1)
            {
              theWalkMethod = &DiffusionProcess::volumeWalk;
            }
          else
            {
              theWalkMethod = &DiffusionProcess::volumeWalkCollide;
            }
        }
      else
        {
          if(WalkProbability == 1)
            {
              theWalkMethod = &DiffusionProcess::surfaceWalk;
            }
          else
            {
              theWalkMethod = &DiffusionProcess::surfaceWalkCollide;
            }
        }
    }
  virtual void printParameters()
    {
      for(std::vector<Species*>::const_iterator i(theDiffusionSpecies.begin());
          i != theDiffusionSpecies.end(); ++i)
        {
          std::cout << getIDString(*i) << " ";
        }
      std::cout << ":" << std::endl << "  Diffusion interval=" <<
        theStepInterval << ", D=" << D << ", Walk probability (P/rho)=" <<
        WalkProbability << std::endl;
    }
  virtual void fire()
    {
      (this->*theWalkMethod)();
      theTime += theStepInterval;
      thePriorityQueue->moveTop();
    }
  void volumeWalk() const
    {
      for(std::vector<Species*>::const_iterator i(theDiffusionSpecies.begin());
          i != theDiffusionSpecies.end(); ++i)
        {
          (*i)->volumeWalk();
        }
    }
  void volumeWalkCollide() const
    {
      for(std::vector<Species*>::const_iterator i(theDiffusionSpecies.begin());
          i != theDiffusionSpecies.end(); ++i)
        {
          (*i)->volumeWalkCollide();
        }
    }
  void surfaceWalk() const
    {
      for(std::vector<Species*>::const_iterator i(theDiffusionSpecies.begin());
          i != theDiffusionSpecies.end(); ++i)
        {
          (*i)->surfaceWalk();
        }
    }
  void surfaceWalkCollide() const
    {
      for(std::vector<Species*>::const_iterator i(theDiffusionSpecies.begin());
          i != theDiffusionSpecies.end(); ++i)
        {
          (*i)->surfaceWalkCollide();
        }
    }
  virtual void initializeLastOnce()
    {
      for(std::vector<Species*>::const_iterator i(theDiffusionSpecies.begin());
          i != theDiffusionSpecies.end(); ++i)
        {
          (*i)->addInterruptedProcess(this);
        }
    }
  virtual void addSubstrateInterrupt(Species* aSpecies, Voxel* aMolecule)
    {
      if(theStepInterval == libecs::INF)
        {
          theStepInterval = theDiffusionSpecies[0]->getDiffusionInterval();
          theTime = theSpatiocyteStepper->getCurrentTime() + theStepInterval; 
          thePriorityQueue->move(theQueueID);
        }
    }
  virtual void removeSubstrateInterrupt(Species* aSpecies, Voxel* aMolecule)
    {
      if(theStepInterval != libecs::INF)
        {
          for(std::vector<Species*>::const_iterator
              i(theDiffusionSpecies.begin());
              i != theDiffusionSpecies.end(); ++i)
            {
              if((*i)->size())
                {
                  return;
                }
            }
          theStepInterval = libecs::INF;
          theTime = theStepInterval; 
          thePriorityQueue->move(theQueueID);
        }
    }
protected:
  bool isVolume;
  double D;
  double P;
  double WalkProbability;
  Species* theVacantSpecies;
  WalkMethod theWalkMethod;
  std::vector<Species*> theDiffusionSpecies;
};

#endif /* __DiffusionProcess_hpp */
