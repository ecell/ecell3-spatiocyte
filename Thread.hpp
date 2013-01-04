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

#define ACCESS_ONCE(x) (*(volatile typeof(x) *)&(x))
#define barrier() __asm__ __volatile__("": : :"memory")
#define FLAG_STOP 0
#define FLAG_RUN  1


class Thread
{ 
public: 
  Thread(unsigned anID, unsigned& aThreadsRunning, char& aFlagA, char& aFlagB,
         SpatiocyteStepper& aStepper):
    theID(anID),
    theStepper(aStepper),
    nThreadsRunning(aThreadsRunning),
    flagA(aFlagA),
    flagB(aFlagB)
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
protected:
   virtual void work()
     {
       while(ACCESS_ONCE(flagA) == FLAG_STOP)
         {
           continue;
         }
       theStepper.allocateLattice(theID);
       __sync_fetch_and_add(&nThreadsRunning, 1);
       while(ACCESS_ONCE(flagB) == FLAG_STOP)
         {
           continue;
         }
       theStepper.constructLattice(theID);
       __sync_fetch_and_add(&nThreadsRunning, 1);
       while(ACCESS_ONCE(flagA) == FLAG_STOP)
         {
           continue;
         }
       theStepper.concatenateLattice(theID);
       __sync_fetch_and_add(&nThreadsRunning, 1);
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
};

#endif /* __Thread_hpp */
