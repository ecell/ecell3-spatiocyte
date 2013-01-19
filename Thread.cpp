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
  if(!theID)
    {
      theStepper.runThreads();
    }
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
}

void Thread::initializeLists()
{
  if(!theID)
    {
      theStepper.runThreads();
    }
  theSpecies[0]->initializeLists(theID, theRng, theMols, theTars,
                                 theAdjMols, theAdjTars, theAdjoins, theIDs,
                                 theAdjBoxes);
}

void Thread::walk()
{
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
  if(!theID)
    {
      theStepper.runThreads();
      //theSpecies[0]->walk(theID, 1, 0, theRng, theMols, theTars, theAdjMols, theAdjTars, theRepeatAdjMols, theRepeatAdjTars, theAdjoins, theIDs, theAdjBoxes);
    }
  theSpecies[0]->walk(theID, r, w, theRng, theMols, theTars, theAdjMols, theAdjTars, theAdjAdjMols, theAdjAdjTars, theBorderMols, theBorderTars, theRepeatAdjMols, theRepeatAdjTars, theAdjoins, theIDs, theAdjBoxes);
}

void Thread::work()
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
  initialize();
  __sync_fetch_and_add(&nThreadsRunning, 1);
  while(ACCESS_ONCE(flagB) == FLAG_STOP)
    {
      continue;
    }
  initializeLists();
  __sync_fetch_and_add(&nThreadsRunning, 1);
  for(;;)
    {
      while(ACCESS_ONCE(flagA) == FLAG_STOP)
        {
          continue;
        }
      __sync_fetch_and_add(&nThreadsRunning, 1);
      walk();
      while(ACCESS_ONCE(flagB) == FLAG_STOP)
        {
          continue;
        }
      __sync_fetch_and_add(&nThreadsRunning, 1);
      walk();
    }
}

