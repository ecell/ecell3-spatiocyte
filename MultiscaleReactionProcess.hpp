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


#ifndef __MultiscaleReactionProcess_hpp
#define __MultiscaleReactionProcess_hpp

#include <sstream>
#include "DiffusionInfluencedReactionProcess.hpp"
#include "SpatiocyteSpecies.hpp"

LIBECS_DM_CLASS(MultiscaleReactionProcess, DiffusionInfluencedReactionProcess)
{ 
public:
  LIBECS_DM_OBJECT(MultiscaleReactionProcess, Process)
    {
      INHERIT_PROPERTIES(DiffusionInfluencedReactionProcess);
    }
  MultiscaleReactionProcess():
    isFinalized(true) {}
  virtual ~MultiscaleReactionProcess() {}
  virtual void initializeThird()
    { 
      theSubstrates.resize(2);
      theProducts.resize(2);
      addedMols.resize(2);
      removedMols.resize(2);
      //This must be in initializeSecond or later since we need to know
      //if a species is multiscale, which is only set by the
      //CompartmentProcess in initializeFirst:
      if(A->getIsMultiscale())
        {
          theMultiscale = A;
          theSubstrates[0] = B;
          theProducts[1] = B;
        }
      else if(B->getIsMultiscale())
        {
          theMultiscale = B;
          theSubstrates[0] = A;
          theProducts[1] = A;
        }
      else
        {
          THROW_EXCEPTION(ValueError, String(
             getPropertyInterface().getClassName()) + " [" + 
              getFullID().asString() + "]: This process must have at least " +
             "one multiscale substrate species.");
        }
      if(theSubstrates[0]->getVacantSpecies() == theMultiscale)
        {
          THROW_EXCEPTION(ValueError, String(
             getPropertyInterface().getClassName()) + " [" + 
              getFullID().asString() + "]: The substrate " + 
              getIDString(theSubstrates[0]) + "'s vacant species is " +
              getIDString(theMultiscale) + " which is a multiscale species. " +
              "This reaction only expects the product's vacant species to be " +
              "a multiscale species. You should probably invert the " +
              "substrate with the product species to reverse the reaction.");
        }
      if(!D)
        {
          THROW_EXCEPTION(ValueError, String(
             getPropertyInterface().getClassName()) + " [" + 
              getFullID().asString() + "]: This process must have two " +
              "products.");
        }
      if(C->getIsMultiscale() && !D->getIsMultiscale())
        {
          theProducts[0] = D;
          theSubstrates[1] = D;
        }
      else if(!C->getIsMultiscale() && D->getIsMultiscale())
        {
          theProducts[0] = C;
          theSubstrates[1] = C;
        }
      else
        {
          THROW_EXCEPTION(ValueError, String(
             getPropertyInterface().getClassName()) + " [" + 
              getFullID().asString() + "]: This process must have at least " +
             "one multiscale product species.");
        }
      //This must be set in
      //initializeThird since it requires vacant species properties
      //set by DiffusionProcess in initializeSecond:

      //If it is a dissociation reaction,
      //theSubstrate diffuses on theMultiscale,
      //theSubstrate unbinds from theMultiscale to become theProduct:
      theMultiscale->setMultiscaleBindIDs(theSubstrates[0]->getID(),
                                            theProducts[0]->getID());
      theMultiscale->setMultiscaleUnbindIDs(theSubstrates[1]->getID(),
                                            theProducts[1]->getID());
      theMultiscale->setDiffusionInfluencedReaction(
            dynamic_cast<DiffusionInfluencedReactionProcess*>(this),
            theSubstrates[0]->getID(), 1); 
    }
  virtual void react(Voxel* aVoxel, const unsigned coord,
                     const unsigned dirA, const unsigned dirB)
    {
      /*
      removedMols[dirA].push_back(aVoxel->coord);
      aVoxel->id = theProducts[dirA]->getID();
      for(unsigned i(0); i != removedMols[dirB].size(); ++i)
        {
          if(removedMols[dirB][i] == coord)
            {
              return;
            }
        }
      addedMols[dirA].push_back(coord);
      */
      theSubstrates[dirA]->softRemoveMolecule(aVoxel);
      theProducts[dirA]->addMolecule(aVoxel);
    }
  virtual void finalizeReaction()
    {
      /*
      theSubstrates[0]->updateMoleculeList(addedMols[1]);
      theSubstrates[1]->updateMoleculeList(addedMols[0]);
      addedMols[0].resize(0);
      addedMols[1].resize(0);
      removedMols[0].resize(0);
      removedMols[1].resize(0);
      for(unsigned k(0); k != 2; ++k)
        {
          //theSubstrates[k]->updateMoleculeList();
          for(unsigned i(0); i != theSubstrates[k]->size(); ++i)
            {
              for(unsigned j(0); j != theSubstrates[k]->size(); ++j)
                {
                  if(i != j && theSubstrates[k]->getMolecule(i) == 
                     theSubstrates[k]->getMolecule(j))
                    {
                      std::cout << " after update error in:" << getFullID().asString() << " sp:" << theSubstrates[k]->getIDString() << " time:" << getStepper()->getCurrentTime() << " i:" << i << " j:" << j << " size:" << theSubstrates[k]->size() << std::endl;
                    }
                }
            }
        }
        */
      DiffusionInfluencedReactionProcess::finalizeReaction();
    }
protected:
  std::vector<Species*> theSubstrates;
  std::vector<Species*> theProducts;
  std::vector<std::vector<unsigned> > addedMols;
  std::vector<std::vector<unsigned> > removedMols;
  Species* theMultiscale;
  double time;
  bool isFinalized;
};

#endif /* __MultiscaleReactionProcess_hpp */

