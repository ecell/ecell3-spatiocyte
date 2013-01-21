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
#include "Species.hpp"

#define ACCESS_ONCE(x) (*(volatile typeof(x) *)&(x))
#define barrier() __asm__ __volatile__("": : :"memory")
#define FLAG_STOP 0
#define FLAG_RUN  1

class Thread
{ 
public: 
  Thread(unsigned anID, unsigned aThreadSize, unsigned& aThreadsRunning,
         char& aFlagA, char& aFlagB, std::vector<Species*>& aSpecies,
         SpatiocyteStepper& aStepper):
    theID(anID),
    theThreadSize(aThreadSize),
    theStepper(aStepper),
    nThreadsRunning(aThreadsRunning),
    flagA(aFlagA),
    flagB(aFlagB),
    theSpecies(aSpecies),
    isToggled(false),
    isRunA(true)
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
  unsigned getBorderMolsSize(unsigned r)
    {
      unsigned aSize(0);
      for(unsigned i(0); i != theBoxSize; ++i)
        {
          aSize += theBorderMols[r][i].size();
        }
      return aSize;
    }
  unsigned getBorderTarsSize(unsigned r)
    {
      unsigned aSize(0);
      for(unsigned i(0); i != theBoxSize; ++i)
        {
          aSize += theBorderTars[r][i].size();
        }
      return aSize;
    }
  unsigned getAdjMolsSize(unsigned r)
    {
      unsigned aSize(0);
      for(unsigned i(0); i != theBoxSize; ++i)
        {
          aSize += theAdjMols[r][i].size();
        }
      return aSize;
    }
  unsigned getAdjTarsSize(unsigned r)
    {
      unsigned aSize(0);
      for(unsigned i(0); i != theBoxSize; ++i)
        {
          aSize += theAdjTars[r][i].size();
        }
      return aSize;
    }
  unsigned getAdjAdjMolsSize(unsigned r)
    {
      unsigned aSize(0);
      for(unsigned i(0); i != theBoxSize; ++i)
        {
          aSize += theAdjAdjMols[r][i].size();
        }
      return aSize;
    }
  unsigned getAdjAdjTarsSize(unsigned r)
    {
      unsigned aSize(0);
      for(unsigned i(0); i != theBoxSize; ++i)
        {
          aSize += theAdjAdjTars[r][i].size();
        }
      return aSize;
    }
  unsigned getRepeatAdjMolsSize()
    {
      unsigned aSize(0);
      for(unsigned i(0); i != theBoxSize; ++i)
        {
          aSize += theRepeatAdjMols[i].size();
        }
      return aSize;
    }
  unsigned getRepeatAdjTarsSize()
    {
      unsigned aSize(0);
      for(unsigned i(0); i != theBoxSize; ++i)
        {
          aSize += theRepeatAdjTars[i].size();
        }
      return aSize;
    }
  std::vector<unsigned>& getBorderMols(unsigned aBox, unsigned r)
    {
      return theBorderMols[r][aBox];
    }
  std::vector<unsigned>& getBorderTars(unsigned aBox, unsigned r)
    {
      return theBorderTars[r][aBox];
    }
  std::vector<unsigned>& getAdjMols(unsigned aBox, unsigned r)
    {
      return theAdjMols[r][aBox];
    }
  std::vector<unsigned>& getAdjTars(unsigned aBox, unsigned r)
    {
      return theAdjTars[r][aBox];
    }
  std::vector<unsigned>& getAdjAdjMols(unsigned aBox, unsigned r)
    {
      return theAdjAdjMols[r][aBox];
    }
  std::vector<unsigned>& getAdjAdjTars(unsigned aBox, unsigned r)
    {
      return theAdjAdjTars[r][aBox];
    }
  void pushAdj(unsigned aBox, unsigned rw, unsigned aMol, unsigned aTar)
    {
      theAdjMols[rw][aBox].push_back(aMol);
      theAdjTars[rw][aBox].push_back(aTar);
    }
  void setMolID(unsigned aMol, unsigned short anID)
    {
      theIDs[aMol] = anID;
    }
  void updateMols(std::vector<unsigned>& aMols);
  void runChildren();
  void waitChildren();
  void waitParent();
  void initialize();
  void initializeLists();
  void walk();
protected:
  void work();
  static void* enter(void* arg)
    {
      ((Thread*)arg)->work();
      return NULL;
    }
  pthread_t theThreadID;
  const unsigned theID;
  const unsigned theThreadSize;
  SpatiocyteStepper& theStepper;
  unsigned& nThreadsRunning;
  char& flagA;
  char& flagB;
  std::ofstream out;
  std::vector<Species*>& theSpecies;
  std::vector<unsigned> theTars;
  std::vector<unsigned> theMols;
  std::vector<unsigned> theAdjoins;
  std::vector<unsigned short> theIDs;
  std::vector<unsigned> theAdjBoxes;
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
  bool isToggled;
  unsigned startBox;
  unsigned endBox;
  bool isRunA;
};

#endif /* __Thread_hpp */
