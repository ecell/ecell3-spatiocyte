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

#include "Species.hpp"

//all writing -> local
//all reading -> shared
void Species::updateBoxMols(const unsigned currBox, const unsigned r,
                            std::vector<unsigned>& aMols,
                            std::vector<unsigned>& aTars,
                            const std::vector<unsigned>& anAdjBoxes)
{
  for(unsigned i(0); i != anAdjBoxes.size(); ++i)
    {
      const unsigned adjBox(anAdjBoxes[i]);
      //reading border mols, so shared:
      std::vector<unsigned>& borderMols(
                  theThreads[adjBox]->getBorderMols(currBox, r));
      std::vector<unsigned>& borderTars(
                  theThreads[adjBox]->getBorderTars(currBox, r));
      for(unsigned j(0); j < borderMols.size(); ++j)
        {
          aMols.push_back(borderMols[j]);
          aTars.push_back(borderTars[j]);
        }
      borderMols.resize(0);
      borderTars.resize(0);
    }
}

void Species::walkMols(std::vector<unsigned>& aMols,
                       const std::vector<unsigned>& aTars,
                       std::vector<unsigned short>& anIDs)
{
  for(unsigned i(0); i < aMols.size(); ++i)
    {
      const unsigned aTar(aTars[i]);
      const unsigned aTarMol(aTar%theBoxMaxSize);
      if(anIDs[aTarMol] == theVacantID)
        {
          if(theWalkProbability == 1 ||
             gsl_rng_uniform(theRng) < theWalkProbability)
            {
              anIDs[aTarMol] = theID;
              anIDs[aMols[i]] = theVacantID;
              aMols[i] = aTarMol;
            }
        }
    }
}

void Species::updateAdjMols(const unsigned currBox, const unsigned r,
                            std::vector<std::vector<unsigned> >& aRepeatAdjMols,
                            std::vector<std::vector<unsigned> >& aRepeatAdjTars,
                            /*
                            std::vector<std::vector<unsigned> >& anAdjMols,
                            std::vector<std::vector<unsigned> >& anAdjTars,
                            */
                            const std::vector<unsigned>& anAdjBoxes)
{
  for(unsigned i(0); i != anAdjBoxes.size(); ++i)
    {
      const unsigned adjBox(anAdjBoxes[i]);
      std::vector<unsigned>& repeatAdjMols(aRepeatAdjMols[adjBox]);
      std::vector<unsigned>& repeatAdjTars(aRepeatAdjTars[adjBox]);
      /*
      std::vector<unsigned>& adjMols(anAdjMols[adjBox]);
      std::vector<unsigned>& adjTars(anAdjTars[adjBox]);
      */
      std::vector<unsigned>& adjMols(theThreads[adjBox
                                     ]->getAdjMols(currBox, r));
      std::vector<unsigned>& adjTars(theThreads[adjBox
                                     ]->getAdjTars(currBox, r));
      for(unsigned j(0); j != repeatAdjMols.size(); ++j)
        {
          adjMols.push_back(repeatAdjMols[j]);
          adjTars.push_back(repeatAdjTars[j]);
        }
      repeatAdjMols.resize(0);
      repeatAdjTars.resize(0);
    }
}

void Species::updateAdjAdjMols(const unsigned currBox, const unsigned r)
{
  //for(unsigned i(0); i != anAdjAdjBoxes.size(); ++i)
  for(unsigned i(0); i != theBoxSize; ++i)
    {
      //const unsigned adjAdjBox(anAdjAdjBoxes[i]);
      const unsigned adjAdjBox(i);
      std::vector<unsigned>& adjAdjMols(theThreads[adjAdjBox
                                        ]->getAdjAdjMols(currBox, r));
      std::vector<unsigned>& adjAdjTars(theThreads[adjAdjBox
                                        ]->getAdjAdjTars(currBox, r));
      for(unsigned j(0); j != adjAdjMols.size(); ++j)
        {
          const unsigned aMol(adjAdjMols[j]);
          const unsigned aBox(aMol/theBoxMaxSize);
          theThreads[aBox]->pushAdj(currBox, r, aMol-theBoxMaxSize*aBox,
                                   adjAdjTars[j]);
        }
      adjAdjMols.resize(0);
      adjAdjTars.resize(0);
    }
}

void Species::setTars(const unsigned currBox,
                      const unsigned w,
                      std::vector<unsigned>& aMols,
                      std::vector<unsigned>& aTars,
                      std::vector<std::vector<unsigned> >& anAdjMols,
                      std::vector<std::vector<unsigned> >& anAdjTars,
                      const std::vector<unsigned>& anAdjoins,
                      RandomLib::Random& aRng)
{
  aTars.resize(0);
  for(unsigned i(0); i < aMols.size(); ++i)
    {
      unsigned& aMol(aMols[i]);
      const unsigned aTar(anAdjoins[aMol*theAdjoinSize+aRng.IntegerC(11)]);
      if(aTar/theBoxMaxSize == currBox) 
        {
          aTars.push_back(aTar);
        }
      else
        {
          //theThreads[aTar/theBoxMaxSize]->pushAdj(currBox, w, aMol, aTar);
          //theThreads[aTar/theBoxMaxSize/(theBoxSize/theThreads.size())]->pushAdj(currBox, w, aMol, aTar);
          anAdjMols[aTar/theBoxMaxSize].push_back(aMol);
          anAdjTars[aTar/theBoxMaxSize].push_back(aTar);
          aMol = aMols.back();
          aMols.pop_back();
          --i;
        }
    }
}
void Species::setAdjTars(const unsigned currBox, const unsigned r,
                std::vector<std::vector<unsigned> >& aBorderMols,
                std::vector<std::vector<unsigned> >& aBorderTars,
                std::vector<std::vector<unsigned> >& anAdjAdjMols,
                std::vector<std::vector<unsigned> >& anAdjAdjTars,
                std::vector<std::vector<unsigned> >& aRepeatAdjMols,
                std::vector<std::vector<unsigned> >& aRepeatAdjTars,
                const std::vector<unsigned>& anAdjBoxes,
                const std::vector<unsigned>& anAdjoins,
                RandomLib::Random& aRng)

{
  for(unsigned i(0); i != anAdjBoxes.size(); ++i)
    {
      const unsigned adjBox(anAdjBoxes[i]);
      //reading adjMols, so get it from the thread:
      std::vector<unsigned>& adjMols(theThreads[adjBox
                                     ]->getAdjMols(currBox, r));
      for(unsigned j(0); j < adjMols.size(); ++j)
        {
          const unsigned aMol(adjMols[j]);
          const unsigned aTar(anAdjoins[aMol*theAdjoinSize+aRng.IntegerC(11)]);
          const unsigned aBox(aTar/theBoxMaxSize);
          if(aBox == currBox) 
            {
              aRepeatAdjMols[adjBox].push_back(aMol);
              aRepeatAdjTars[adjBox].push_back(aTar);
            }
          else if(aBox == adjBox)
            {
              aBorderMols[adjBox].push_back(aMol);
              aBorderTars[adjBox].push_back(aTar);
            }
          else
            {
              anAdjAdjMols[aBox].push_back(theBoxMaxSize*adjBox+aMol);
              anAdjAdjTars[aBox].push_back(aTar);
            }
        }
      adjMols.resize(0);
    }
}

void Species::walkAdjMols(const unsigned currBox, const unsigned r,
                          std::vector<unsigned>& aMols,
                          std::vector<unsigned short>& anIDs,
                          const std::vector<unsigned>& anAdjBoxes)
{
  for(unsigned i(0); i != anAdjBoxes.size(); ++i)
    {
      const unsigned aBox(anAdjBoxes[i]);
      std::vector<unsigned>& adjMols(theThreads[aBox]->getAdjMols(currBox, r));
      std::vector<unsigned>& adjTars(theThreads[aBox]->getAdjTars(currBox, r));
      for(unsigned j(0); j < adjMols.size(); ++j)
        {
          const unsigned aTar(adjTars[j]);
          const unsigned aTarMol(aTar%theBoxMaxSize);
          if(anIDs[aTarMol] == theVacantID)
            {
              anIDs[aTarMol] = theID;
              theThreads[aBox]->setMolID(adjMols[j], theVacantID);
              //theIDs[aBox][adjMols[j]] = theVacantID;
              aMols.push_back(aTarMol);
              adjMols[j] = adjMols.back();
              adjMols.pop_back();
              adjTars[j] = adjTars.back();
              adjTars.pop_back();
              --j;
            }
        }
      adjTars.resize(0);
    }
}


void Species::walk(const unsigned anID, unsigned r, unsigned w,
           RandomLib::Random& aRng,
           std::vector<unsigned>& aMols,
           std::vector<unsigned>& aTars,
           std::vector<std::vector<std::vector<unsigned> > >& anAdjMols,
           std::vector<std::vector<std::vector<unsigned> > >& anAdjTars,
           std::vector<std::vector<std::vector<unsigned> > >& anAdjAdjMols,
           std::vector<std::vector<std::vector<unsigned> > >& anAdjAdjTars,
           std::vector<std::vector<std::vector<unsigned> > >& aBorderMols,
           std::vector<std::vector<std::vector<unsigned> > >& aBorderTars,
           std::vector<std::vector<unsigned> >& aRepeatAdjMols,
           std::vector<std::vector<unsigned> >& aRepeatAdjTars,
           std::vector<unsigned>& anAdjoins,
           std::vector<unsigned short>& anIDs,
           std::vector<unsigned>& anAdjBoxes)
{
  updateBoxMols(anID, r, aMols, aTars, anAdjBoxes);
  walkMols(aMols, aTars, anIDs);
  updateAdjMols(anID, r, aRepeatAdjMols, aRepeatAdjTars, anAdjBoxes);
  updateAdjAdjMols(anID, r); 
  walkAdjMols(anID, r, aMols, anIDs, anAdjBoxes);
  setAdjTars(anID, r, aBorderMols[w], aBorderTars[w], anAdjAdjMols[w],
             anAdjAdjTars[w], aRepeatAdjMols, aRepeatAdjTars, anAdjBoxes,
             anAdjoins, aRng);
  setTars(anID, w, aMols, aTars, anAdjMols[r], anAdjTars[r], anAdjoins, aRng);
  //setTars(theMols[anID], theTars[anID], theAdjMols[w][anID], theAdjTars[w][anID], anID, theAdjoins[anID], aRng);
  //setTars(aMols, aTars, anID, anAdjoins, aRng);
  //setTars(theMols[anID], theTars[anID], anID, theAdjoins[anID], aRng);
  /*
  for(unsigned i(0); i != 6; ++i)
    {
      setTars(theMols[i], theTars[i], anID, theAdjoins[i], aRng);
    }
    */
  /*
  if(!anID)
    {
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
      theStepper->runThreads();
    }
  //for(unsigned i(0); i != theBoxSize; ++i)
  for(unsigned i(anID*2); i != (anID*2)+2; ++i)
    {
      updateBoxMols(theBorderMols[r][i], theBorderTars[r][i], theMols[i],
                    theTars[i], theAdjBoxes[i]);
    }
  for(unsigned i(anID*2); i != (anID*2)+2; ++i)
    {
      walkMols(theMols[i], theTars[i], theIDs[i]);
    }
  for(unsigned i(anID*2); i != (anID*2)+2; ++i)
    {
      updateAdjMols(theRepeatAdjMols[i], theRepeatAdjTars[i],
                    theAdjMols[r][i], theAdjTars[r][i], theAdjBoxes[i]);
    }
  for(unsigned i(anID*2); i != (anID*2)+2; ++i)
    {
      updateAdjAdjMols(theAdjAdjMols[r][i], theAdjAdjTars[r][i],
                    theAdjMols[r][i], theAdjTars[r][i]);
    }
  for(unsigned i(anID*2); i != (anID*2)+2; ++i)
    {
      walkAdjMols(theMols[i], theAdjMols[r][i], theAdjTars[r][i],
                  theIDs[i], theAdjBoxes[i]);
    }
  for(unsigned i(anID*2); i != (anID*2)+2; ++i)
    {
      setAdjTars(theBorderMols[w], theBorderTars[w], theAdjAdjMols[w],
                 theAdjAdjTars[w], theRepeatAdjMols[i],
                 theRepeatAdjTars[i], theAdjMols[r][i], theAdjBoxes[i], i);
    }
  for(unsigned i(anID*2); i != (anID*2)+2; ++i)
    {
      setTars(theMols[i], theTars[i], theAdjMols[w], theAdjTars[w], i,
              theAdjoins[i]);
    }
    */
}

void Species::updateMols()
{
  if(isDiffusiveVacant || isReactiveVacant)
    {
      updateVacantMols();
    }
  else if(isTag)
    {
      updateTagMols();
    }
  if(!theID)
    {
      for(unsigned i(0); i != theThreads.size(); ++i)
        {
          theThreads[i]->updateMols(theMols[i]);
        }
    }
}

