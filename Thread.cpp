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
  //doWork();
  waitChildren();
}

void Thread::initializeLists()
{
  runChildren();
  theSpecies[0]->initializeLists(theID, theRng, theMols, theTars,
                                 theAdjMols, theAdjTars, theAdjoins, theIDs,
                                 theAdjBoxes, theAdjAdjBoxes, theRands);
  waitChildren();
}

void Thread::doWork()
{
  for(unsigned i(0); i != theBoxSize; ++i)
    {
      for(unsigned j(0); j != 20000; ++j)
        {
          theBorderMols[0][i].push_back(theRng.IntegerC(11));
          theBorderMols[1][i].push_back(theRng.IntegerC(11));
          theBorderTars[0][i].push_back(theRng.IntegerC(11));
          theBorderTars[1][i].push_back(theRng.IntegerC(11));
        }
    }
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
  //if(!theID)
    {
  theSpecies[0]->walk(theID, r, w, theRng, theMols, theTars, theAdjMols, theAdjTars, theAdjAdjMols, theAdjAdjTars, theBorderMols, theBorderTars, theRepeatAdjMols, theRepeatAdjTars, theAdjoins, theIDs, theAdjBoxes, theAdjAdjBoxes, theRands);
    }
  waitChildren();
}
void Thread::walk(std::vector<std::vector<std::vector<unsigned> > >& aBorderMols, std::vector<std::vector<std::vector<unsigned> > >& aBorderTars)
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
  //theSpecies[0]->walk(theID, r, w, theRng, theMols, theTars, theBorderMols, theBorderTars, theAdjBoxes);
  //theSpecies[0]->walk(theID, r, w, theRng, theMols, theTars, theAdjMols, theAdjTars, theAdjAdjMols, theAdjAdjTars, aBorderMols, aBorderTars, theRepeatAdjMols, theRepeatAdjTars, theAdjoins, theIDs, theAdjBoxes);
}

void Thread::runChildren()
{
  if(!theID)
    {
      out << "run val:" << nThreadsRunning << std::endl << std::flush;
      __sync_synchronize();
      if(isRunA)
        {
          out << "run flagA" << std::endl << std::flush;
          flagA = FLAG_RUN;
        }
      else
        {
          out << "run flagB" << std::endl << std::flush;
          flagB = FLAG_RUN;
        }
    }
}

void Thread::waitChildren()
{
  if(!theID)
    {
      out << "wait val:" << nThreadsRunning << std::endl << std::flush;
      __sync_synchronize();
      while(ACCESS_ONCE(nThreadsRunning) < theThreadSize-1)
        {
          continue;
        }
      nThreadsRunning = 0;
      __sync_synchronize();
      if(isRunA)
        {
          out << "stop flagA" << std::endl << std::flush;
          flagA = FLAG_STOP;
          isRunA = false;
        }
      else
        {
          out << "stop flagB" << std::endl << std::flush;
          flagB = FLAG_STOP;
          isRunA = true;
        }
    }
}

void Thread::waitParent()
{
  if(isRunA)
    {
      out << "wait flagA:" << nThreadsRunning << std::endl << std::flush;
      __sync_synchronize();
      while(ACCESS_ONCE(flagA) == FLAG_STOP)
        {
          continue;
        }
      isRunA = false;
    }
  else
    {
      out << "wait flagB:" << nThreadsRunning << std::endl << std::flush;
      __sync_synchronize();
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
  out << "val before:" << __sync_fetch_and_add(&nThreadsRunning, 1) << std::endl << std::flush;
  out << "val after:" << nThreadsRunning << std::endl << std::flush;
  waitParent();
  theStepper.concatenateLattice(theID);
  out << "val before:" << __sync_fetch_and_add(&nThreadsRunning, 1) << std::endl << std::flush;
  out << "val after:" << nThreadsRunning << std::endl << std::flush;
  waitParent();
  initialize();
  out << "val before:" << __sync_fetch_and_add(&nThreadsRunning, 1) << std::endl << std::flush;
  out << "val after:" << nThreadsRunning << std::endl << std::flush;
  waitParent();
  initializeLists();
  out << "val before:" << __sync_fetch_and_add(&nThreadsRunning, 1) << std::endl << std::flush;
  out << "val after:" << nThreadsRunning << std::endl << std::flush;
  out << "1" << std::endl << std::flush;
  waitParent();
  out << "2" << std::endl << std::flush;
  /*
  std::vector<std::vector<std::vector<unsigned> > > aBorderMols;
  std::vector<std::vector<std::vector<unsigned> > > aBorderTars;
  aBorderMols.resize(2);
  aBorderTars.resize(2);
  for(unsigned i(0); i != 2; ++i)
    {
      aBorderMols[i].resize(theBoxSize);
      aBorderTars[i].resize(theBoxSize);
    }
  for(unsigned i(0); i != theBoxSize; ++i)
    {
      for(unsigned j(0); j != 20000; ++j)
        {
          aBorderMols[0][i].push_back(theRng.IntegerC(11));
          aBorderMols[1][i].push_back(theRng.IntegerC(11));
          aBorderTars[0][i].push_back(theRng.IntegerC(11));
          aBorderTars[1][i].push_back(theRng.IntegerC(11));
        }
    }
    */
  for(;;)
    {
  out << "3" << std::endl << std::flush;
      //walk(aBorderMols, aBorderTars);
      walk();
  out << "4" << std::endl << std::flush;
      //walk(aBorderMols, aBorderTars);
  out << "val before:" << __sync_fetch_and_add(&nThreadsRunning, 1) << std::endl << std::flush;
  out << "val after:" << nThreadsRunning << std::endl << std::flush;
  out << "5" << std::endl << std::flush;
      waitParent();
    }
}


void Thread::updateMols(std::vector<unsigned>& aMols)
{
  aMols.resize(theMols.size());
  for(unsigned i(0); i != theMols.size(); ++i)
    {
      aMols[i] = theMols[i];
    }
  for(unsigned i(0); i != theBoxSize; ++i)
    {
      for(unsigned j(0); j != theAdjMols[0][i].size(); ++j)
        {
          aMols.push_back(theAdjMols[0][i][j]);
        }
      for(unsigned j(0); j != theAdjMols[1][i].size(); ++j)
        {
          aMols.push_back(theAdjMols[1][i][j]);
        }
    }
}
