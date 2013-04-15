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

#include "LifetimeLogProcess.hpp"

void LifetimeLogProcess::initialize()
{
  if(isInitialized)
    {
      return;
    }
  SpatiocyteProcess::initialize();
  isPriorityQueued = true;
  for(VariableReferenceVector::iterator
      i(theVariableReferenceVector.begin());
      i != theVariableReferenceVector.end(); ++i)
    {
      Species* aSpecies(theSpatiocyteStepper->variable2species(
                                     (*i).getVariable())); 
      if((*i).getCoefficient() == -1)
        {
          theTrackedSpeciesList.push_back(aSpecies);
        }
      else if((*i).getCoefficient() == 1)
        {
          theUntrackedSpeciesList.push_back(aSpecies);
        }
    }
  if(!theTrackedSpeciesList.size())
    {
      THROW_EXCEPTION(ValueError, String(
                      getPropertyInterface().getClassName()) +
                      "[" + getFullID().asString() + 
                      "]: Lifetime logging requires at least one " +
                      "nonHD variable reference with -1 " +
                      "coefficient as the tracked species, " +
                      "but none is given."); 
    }
  if(!theUntrackedSpeciesList.size())
    {
      THROW_EXCEPTION(ValueError, String(
                      getPropertyInterface().getClassName()) +
                      "[" + getFullID().asString() + 
                      "]: Lifetime logging requires at least one " +
                      "nonHD variable reference with 1 " +
                      "coefficient as the untracked species, " +
                      "but none is given."); 
    }
  for(unsigned i(0); i != theTrackedSpeciesList.size(); ++i)
    {
      theTrackedSpeciesList[i]->setIsTagged();
    }
  LogInterval = libecs::INF;
}

void LifetimeLogProcess::initializeFirst()
{
  IteratingLogProcess::initializeFirst();
  isBindingSite.resize(getStepper()->getProcessVector().size(), false);
}

void LifetimeLogProcess::interruptedPre(ReactionProcess* aProcess)
{
  //std::cout << "Pre aProcess:" << aProcess->getFullID().asString() << std::endl;
  if(aProcess->getA())
    {
      if(logTrackedMolecule(aProcess->getID(), aProcess->getA(),
                            aProcess->getMoleculeA()))
        {
          return;
        }
    }
  if(aProcess->getB())
    {
      if(logTrackedMolecule(aProcess->getID(), aProcess->getB(),
                            aProcess->getMoleculeB()))
        {
          return;
        }
    }
}

void LifetimeLogProcess::interruptedPost(ReactionProcess* aProcess)
{
  //std::cout << "Post aProcess:" << aProcess->getFullID().asString() << std::endl;
  if(aProcess->getC())
    {
      if(initTrackedMolecule(aProcess->getC()))
        {
          return;
        }
    }
 if(aProcess->getD())
    {
      if(initTrackedMolecule(aProcess->getD()))
        {
          return;
        }
    }
}

bool LifetimeLogProcess::logTrackedMolecule(const unsigned anID,
                                            Species* aSpecies,
                                            const Voxel* aMolecule)
{
  for(unsigned i(0); i != theTrackedSpeciesList.size(); ++i)
    {
      if(aSpecies == theTrackedSpeciesList[i])
        {
          if(isBindingSite[anID])
            {
              //std::cout << "is binding site" << std::endl;
            }
          const unsigned anIndex(aSpecies->getIndex(aMolecule));
          const Point aPoint(aSpecies->getPoint(anIndex));
          const Point anOrigin(aSpecies->coord2point(
                                         aSpecies->getTag(anIndex).origin));
          //std::cout << "dist:" << distance(aPoint, anOrigin) << std::endl;
          return true;
        }
    }
  return false;
}

bool LifetimeLogProcess::initTrackedMolecule(Species* aSpecies)
{
  for(unsigned i(0); i != theTrackedSpeciesList.size(); ++i)
    {
      if(aSpecies == theTrackedSpeciesList[i])
        {
          const unsigned anIndex(aSpecies->size()-1);
          Tag& aTag(aSpecies->getTag(anIndex));
          aTag.origin = aSpecies->getCoord(anIndex);
          //std::cout << "origin:" << aTag.origin << std::endl;
          return true;
        }
    }
  return false;
}

bool LifetimeLogProcess::isDependentOnPre(const ReactionProcess* aProcess)
{
  const VariableReferenceVector& aVariableReferences(
                                       aProcess->getVariableReferenceVector()); 
  for(unsigned i(0); i != theTrackedSpeciesList.size(); ++i)
    {
      if(isInVariableReferences(aVariableReferences, -1,
                                theTrackedSpeciesList[i]->getVariable()))
        {
          for(unsigned j(0); j != theUntrackedSpeciesList.size(); ++j)
            {
              if(isInVariableReferences(aVariableReferences, 1,
                                    theUntrackedSpeciesList[j]->getVariable()))
                {
                  for(unsigned k(0); k != theTrackedSpeciesList.size(); ++k)
                    {
                      if(isInVariableReferences(aVariableReferences, 1,
                                    theTrackedSpeciesList[k]->getVariable())) 
                        {
                          isBindingSite[aProcess->getID()] = true;
                          return true;
                        }
                    }
                  return true;
                }
            }
        }
    }
  return false;
}

bool LifetimeLogProcess::isInVariableReferences(const VariableReferenceVector&
                                                aVariableReferences,
                                                const int aCoefficient,
                                                const Variable* aVariable) const
{
  for(VariableReferenceVector::const_iterator
      i(aVariableReferences.begin()); i != aVariableReferences.end(); ++i)
    {
      //If the both coefficients have the same sign:
      if((*i).getCoefficient()*aCoefficient > 0 &&
         (*i).getVariable() == aVariable)
        {
          return true;
        }
    }
  return false;
}

bool LifetimeLogProcess::isDependentOnPost(const ReactionProcess* aProcess)
{
  const VariableReferenceVector& aVariableReferences(
                                       aProcess->getVariableReferenceVector()); 
  for(unsigned i(0); i != theTrackedSpeciesList.size(); ++i)
    {
      if(isInVariableReferences(aVariableReferences, 1,
                                theTrackedSpeciesList[i]->getVariable()))
        {
          for(unsigned j(0); j != theTrackedSpeciesList.size(); ++j)
            {
              if(isInVariableReferences(aVariableReferences, -1,
                                    theTrackedSpeciesList[j]->getVariable()))
                {
                  return false;
                }
            }
          return true;
        }
    }
  return false;
}


LIBECS_DM_INIT(LifetimeLogProcess, Process); 
