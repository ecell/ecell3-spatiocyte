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

#include "Thread.hpp"

void Thread::initialize()
{
  runChildren();
  theBoxSize = theStepper.getBoxSize();
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
  waitChildren();
}

void Thread::initializeLists()
{
  runChildren();
  theSpecies[0]->initializeLists(theID, theRng, theMols, theTars,
                                 theAdjMols, theAdjTars, theAdjoins, theIDs,
                                 theAdjBoxes);
  waitChildren();
}

void Thread::walk()
{
  runChildren();
  unsigned r(0);
  unsigned w(1);
  if(isToggled)
    {
      r = 1;
      w = 0;
      isToggled = false;
    }
  else
    {
      isToggled = true;
    } 
  theSpecies[0]->walk(theID, r, w, theRng, theMols, theTars, theAdjMols,
                      theAdjTars, theAdjAdjMols, theAdjAdjTars, theBorderMols,
                      theBorderTars, theRepeatAdjMols, theRepeatAdjTars,
                      theAdjoins, theIDs, theAdjBoxes);
  waitChildren();
}

void Thread::runChildren()
{
  if(!theID)
    {
      if(isRunA)
        {
          flagA = FLAG_RUN;
        }
      else
        {
          flagB = FLAG_RUN;
        }
    }
}

void Thread::waitChildren()
{
  if(!theID)
    {
      barrier();
      while(ACCESS_ONCE(nThreadsRunning) < theThreadSize-1)
        {
          continue;
        }
      nThreadsRunning = 0;
      if(isRunA)
        {
          flagA = FLAG_STOP;
          isRunA = false;
        }
      else
        {
          flagB = FLAG_STOP;
          isRunA = true;
        }
    }
}

void Thread::waitParent()
{
  if(isRunA)
    {
      while(ACCESS_ONCE(flagA) == FLAG_STOP)
        {
          continue;
        }
      isRunA = false;
    }
  else
    {
      while(ACCESS_ONCE(flagB) == FLAG_STOP)
        {
          continue;
        }
      isRunA = true;
    }
}

void Thread::work()
{
  waitParent();
  theStepper.constructLattice(theID);
  __sync_fetch_and_add(&nThreadsRunning, 1);
  waitParent();
  theStepper.concatenateLattice(theID);
  __sync_fetch_and_add(&nThreadsRunning, 1);
  waitParent();
  initialize();
  __sync_fetch_and_add(&nThreadsRunning, 1);
  waitParent();
  initializeLists();
  __sync_fetch_and_add(&nThreadsRunning, 1);
  for(;;)
    {
      waitParent();
      walk();
      __sync_fetch_and_add(&nThreadsRunning, 1);
    }
}

