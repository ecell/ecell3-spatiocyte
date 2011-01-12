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


#ifndef __IteratingLogProcess_hpp
#define __IteratingLogProcess_hpp

#include <fstream> //provides ofstream
#include <math.h>
#include "SpatiocyteProcess.hpp"
#include "SpatiocyteSpecies.hpp"

LIBECS_DM_CLASS(IteratingLogProcess, SpatiocyteProcess)
{ 
public:
  LIBECS_DM_OBJECT(IteratingLogProcess, Process)
    {
      INHERIT_PROPERTIES(Process);
      PROPERTYSLOT_SET_GET(Real, LogDuration);
      PROPERTYSLOT_SET_GET(Real, LogInterval);
      PROPERTYSLOT_SET_GET(Integer, Iterations);
      PROPERTYSLOT_SET_GET(Integer, SaveInterval);
      PROPERTYSLOT_SET_GET(Integer, Centered);
      PROPERTYSLOT_SET_GET(Integer, InContact);
      PROPERTYSLOT_SET_GET(Integer, Survival);
      PROPERTYSLOT_SET_GET(Integer, RebindTime);
      PROPERTYSLOT_SET_GET(Integer, Displacement);
      PROPERTYSLOT_SET_GET(Integer, Diffusion);
      PROPERTYSLOT_SET_GET(String, FileName);
    }
  SIMPLE_SET_GET_METHOD(Real, LogDuration);
  SIMPLE_SET_GET_METHOD(Real, LogInterval);
  SIMPLE_SET_GET_METHOD(Integer, Iterations);
  SIMPLE_SET_GET_METHOD(Integer, SaveInterval);
  SIMPLE_SET_GET_METHOD(Integer, Centered);
  SIMPLE_SET_GET_METHOD(Integer, InContact);
  SIMPLE_SET_GET_METHOD(Integer, Survival);
  SIMPLE_SET_GET_METHOD(Integer, RebindTime);
  SIMPLE_SET_GET_METHOD(Integer, Displacement);
  SIMPLE_SET_GET_METHOD(Integer, Diffusion);
  SIMPLE_SET_GET_METHOD(String, FileName);
  IteratingLogProcess():
    SpatiocyteProcess(false),
    Centered(0),
    Diffusion(0),
    Displacement(0),
    InContact(0),
    RebindTime(0),
    SaveInterval(0),
    Survival(0),
    LogInterval(0) {}
  virtual ~IteratingLogProcess() {}
  virtual void initializeSecond()
    {
      SpatiocyteProcess::initializeSecond(); 
      for(unsigned int i(0); i != theProcessSpecies.size(); ++i)
        {
          if(InContact)
            {
              theProcessSpecies[i]->setIsInContact();
            }
          if(Centered)
            {
              theProcessSpecies[i]->setIsCentered();
            }
        }
      theLogCnt = 0;
    }
  virtual void initializeFourth()
    {
      for(std::vector<Species*>::const_iterator i(theProcessSpecies.begin());
          i != theProcessSpecies.end(); ++i)
        {
          if((*i)->getDiffusionInterval() < theStepInterval)
            {
              theStepInterval = (*i)->getDiffusionInterval();
            }
          Species* reactantPair((*i)->getDiffusionInfluencedReactantPair());
          if(reactantPair != NULL && 
             reactantPair->getDiffusionInterval() < theStepInterval)
            {
              theStepInterval = reactantPair->getDiffusionInterval();
            }
        }
      if(LogInterval > 0)
        {
          theStepInterval = LogInterval;
        }
      else
        {
          LogInterval = theStepInterval;
        }
      theTime = LogInterval;
      thePriorityQueue->move(theQueueID);
    }
  virtual void initializeLastOnce()
    {
      theLogFile.open(FileName.c_str(), std::ios::trunc);
      theTotalIterations = Iterations;
      theLogValues.resize(theProcessSpecies.size());
      for(unsigned int i(0); i != theProcessSpecies.size(); ++i)
        {
          if(RebindTime)
            {
              theLogValues[i].resize(Iterations);
            }
          else
            {
              theLogValues[i].resize((int)ceil(LogDuration/theStepInterval));
            }
          for(unsigned int j(0); j != theLogValues[i].size(); ++j)
            {
              theLogValues[i][j] = 0;
            }
        }
    }
  virtual void fire()
    {
      if(Iterations == 0)
        {
          std::cout << "Saving data in: " << FileName.c_str() << std::endl;
          double aTime(LogInterval);
          for(unsigned int i(0); i != theLogValues[0].size(); ++i)
            {
              theLogFile << std::setprecision(15) << aTime;
              for(unsigned int j(0); j != theProcessSpecies.size(); ++j)
                {
                  if(RebindTime)
                    {
                      theLogFile << "," << std::setprecision(15) << theLogValues[j][i];
                    }
                  else
                    {
                      theLogFile << "," << std::setprecision(15) << theLogValues[j][i]/
                        theTotalIterations;
                    }
                }
              theLogFile << std::endl;
              aTime += LogInterval;
            }
          theLogFile.close();
          theStepInterval = libecs::INF;
        }
      else if(theTime >= LogDuration)
        {
          theStepInterval = LogInterval;
          theSpatiocyteStepper->reset(Iterations);
          --Iterations;
          std::cout << "Iterations left:" << Iterations << std::endl;
          if(SaveInterval > 0 && 
             Iterations%(int)rint(theTotalIterations/SaveInterval) == 0)
            {
              std::string aFileName(FileName.c_str());
              aFileName = aFileName + ".back";
              std::cout << "Saving temporary backup data in: " << aFileName << std::endl;
              std::ofstream aFile;
              aFile.open(aFileName.c_str(), std::ios::trunc);
              double aTime(LogInterval);
              int completedIterations(theTotalIterations-Iterations);
              unsigned int aSize(theLogValues[0].size());
              if(RebindTime)
                {
                  aSize = completedIterations; 
                }
              for(unsigned int i(0); i != aSize; ++i)
                {
                  aFile << std::setprecision(15) << aTime;
                  for(unsigned int j(0); j != theProcessSpecies.size(); ++j)
                    {
                      if(RebindTime)
                        {
                          aFile << "," << std::setprecision(15) <<
                            theLogValues[j][i];
                        }
                      else
                        {
                          aFile << "," << std::setprecision(15) <<
                            theLogValues[j][i]/completedIterations;
                        }
                    }
                  aFile << std::endl;
                  aTime += theStepInterval;
                }
              aFile.close();
            }
          return;
        }
      else
        {
          theSurvivalCnt = 0;
          for(unsigned int i(0); i != theProcessSpecies.size(); ++i)
            {
              if(RebindTime)
                {
                  if(theProcessSpecies[i]->getVariable()->getValue())
                    {
                      ++theSurvivalCnt;
                    }
                  if(!theLogValues[i][theTotalIterations-Iterations])
                    {
                      if(!theProcessSpecies[i]->getVariable()->getValue())
                        {
                          theLogValues[i][theTotalIterations-Iterations] =
                            theTime;
                        }
                    }
                }
              else if(Survival)
                {
                  if(theProcessSpecies[i]->getVariable()->getValue())
                    {
                      ++theSurvivalCnt;
                    }
                  theLogValues[i][theLogCnt] += 
                    theProcessSpecies[i]->getVariable()->getValue()/
                    theProcessSpecies[i]->getInitMoleculeSize();
                }
              else if(Displacement)
                {
                  theLogValues[i][theLogCnt] += 
                    theProcessSpecies[i]->getMeanSquaredDisplacement();
                }
              else if(Diffusion)
                {
                  if(theProcessSpecies[i]->getIsVolume())
                    {
                      theLogValues[i][theLogCnt] += 
                        theProcessSpecies[i]->getMeanSquaredDisplacement()/
                        (6*theTime);
                    }
                  else
                    {
                      theLogValues[i][theLogCnt] += 
                        theProcessSpecies[i]->getMeanSquaredDisplacement()/
                         (4*theTime);
                    }
                }
              //By default log the values:
              else
                {
                  theLogValues[i][theLogCnt] += 
                    theProcessSpecies[i]->getVariable()->getValue();
                }
            }
          //If all survival species are dead, go on to the next iteration:
          if((Survival || RebindTime) && !theSurvivalCnt)
            {
              theStepInterval = LogDuration-theTime; 
            }
          ++theLogCnt;
        }
      theTime += theStepInterval;
      thePriorityQueue->moveTop();
    }
protected:
  int Centered;
  int Diffusion;
  int Displacement;
  int InContact;
  int Iterations;
  int RebindTime;
  int SaveInterval;
  int Survival;
  int theLogCnt;
  int theSurvivalCnt;
  int theTotalIterations;
  double LogDuration;
  double LogInterval;
  String FileName;
  std::ofstream theLogFile;
  std::vector<std::vector<double> > theLogValues;
};

#endif /* __IteratingLogProcess_hpp */
