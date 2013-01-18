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


#ifndef __Thread_hpp
#define __Thread_hpp

#include <fstream>
#include <sstream>
#include <pthread.h>
#include "SpatiocyteStepper.hpp"
#include "SpatiocyteSpecies.hpp"

#define ACCESS_ONCE(x) (*(volatile typeof(x) *)&(x))
#define barrier() __asm__ __volatile__("": : :"memory")
#define FLAG_STOP 0
#define FLAG_RUN  1

class Thread
{ 
public: 
  Thread(unsigned anID, unsigned& aThreadsRunning, char& aFlagA, char& aFlagB,
         std::vector<Species*>& aSpecies, SpatiocyteStepper& aStepper):
    theID(anID),
    theStepper(aStepper),
    nThreadsRunning(aThreadsRunning),
    flagA(aFlagA),
    flagB(aFlagB),
    theSpecies(aSpecies)
  {
    std::ostringstream fileName;
    fileName << "thread" << anID << ".out" << std::ends;
    out.open(fileName.str().c_str());
    __sync_fetch_and_add(&nThreadsRunning, 1);
  }
  virtual ~Thread()
    {
      out.close();
    }
  void create()
    {
      pthread_create(&theThreadID, NULL, enter, this);
    }
  void wait()
    {
      pthread_join(theThreadID, NULL);
    }
  void initializeLists()
    {
      theBoxSize = theStepper.getBoxSize();
      if(!theID)
        {
          theStepper.runThreads();
        }
      theRng.Reseed();
      theAdjMols.resize(2);
      theAdjTars.resize(2);
      theAdjAdjMols.resize(2);
      theAdjAdjTars.resize(2);
      theBorderMols.resize(2);
      theBorderTars.resize(2);
      theRepeatAdjMols.resize(theBoxSize);
      theRepeatAdjTars.resize(theBoxSize);
      for(unsigned i(0); i != 2; ++i)
        {
          theAdjMols[i].resize(theBoxSize);
          theAdjTars[i].resize(theBoxSize);
          theAdjAdjMols[i].resize(theBoxSize);
          theAdjAdjTars[i].resize(theBoxSize);
          theBorderMols[i].resize(theBoxSize);
          theBorderTars[i].resize(theBoxSize);
        }
      theSpecies[0]->initializeLists(theID, theRng, theMols, theTars,
                                     theAdjMols, theAdjTars, theAdjoins);
    }
  void walk()
    {
      if(!theID)
        {
          theStepper.runThreads();
          //setTars(theMols, theTars, theID, theAdjoins, theRng);
          //theSpecies[0]->walk(theID, 1, 0, theRng, theMols, theTars, theAdjMols, theAdjTars, theAdjoins);
        }
      theSpecies[0]->walk(theID, 1, 0, theRng, theMols, theTars, theAdjMols, theAdjTars, theAdjoins);
    }
protected:
  static void setTars(std::vector<unsigned>& aMols,
                      std::vector<unsigned>& aTars,
                      const unsigned currBox,
                      const std::vector<unsigned>& anAdjoins,
                      RandomLib::Random& aRng)
    {
      aTars.resize(aMols.size());
      for(unsigned i(0); i < aMols.size(); ++i)
        {
          unsigned& aMol(aMols[i]);
          const unsigned aTar(anAdjoins[aMol*12+aRng.IntegerC(11)]);
          aTars[i] = aTar;
        }
    }
  virtual void work()
    {
      while(ACCESS_ONCE(flagA) == FLAG_STOP)
        {
          continue;
        }
      theStepper.constructLattice(theID);
      __sync_fetch_and_add(&nThreadsRunning, 1);
      while(ACCESS_ONCE(flagB) == FLAG_STOP)
        {
          continue;
        }
      theStepper.concatenateLattice(theID);
      __sync_fetch_and_add(&nThreadsRunning, 1);
      while(ACCESS_ONCE(flagA) == FLAG_STOP)
        {
          continue;
        }
      initializeLists();
      __sync_fetch_and_add(&nThreadsRunning, 1);
      for(;;)
        {
          while(ACCESS_ONCE(flagB) == FLAG_STOP)
            {
              continue;
            }
          __sync_fetch_and_add(&nThreadsRunning, 1);
          walk();
          /*
          for(unsigned i(0); i != 800; ++i)
            {
              aTars.resize(0);
              for(unsigned j(0); j != 100; ++j)
                {
                  aTars.push_back(j);
                }
            }
            */
          while(ACCESS_ONCE(flagA) == FLAG_STOP)
            {
              continue;
            }
          __sync_fetch_and_add(&nThreadsRunning, 1);
          walk();
          /*
          for(unsigned i(0); i != 800; ++i)
            {
              aTars.resize(0);
              for(unsigned j(0); j != 100; ++j)
                {
                  aTars.push_back(j);
                }
            }
            */
        }
    }
private: 
  static void* enter(void* arg)
    {
      ((Thread*)arg)->work();
      return NULL;
    }
  pthread_t theThreadID;
  unsigned theID;
  SpatiocyteStepper& theStepper;
  unsigned& nThreadsRunning;
  char& flagA;
  char& flagB;
  std::ofstream out;
  std::vector<Species*>& theSpecies;
  std::vector<unsigned> theTars;
  std::vector<unsigned> theMols;
  std::vector<unsigned> theAdjoins;
  RandomLib::Random theRng;
  unsigned theBoxSize;
  std::vector<std::vector<std::vector<unsigned> > > theAdjMols;
  std::vector<std::vector<std::vector<unsigned> > > theAdjTars;
  std::vector<std::vector<std::vector<unsigned> > > theAdjAdjMols;
  std::vector<std::vector<std::vector<unsigned> > > theAdjAdjTars;
  std::vector<std::vector<std::vector<unsigned> > > theBorderMols;
  std::vector<std::vector<std::vector<unsigned> > > theBorderTars;
  std::vector<std::vector<unsigned> > theRepeatAdjMols;
  std::vector<std::vector<unsigned> > theRepeatAdjTars;
};

#endif /* __Thread_hpp */
