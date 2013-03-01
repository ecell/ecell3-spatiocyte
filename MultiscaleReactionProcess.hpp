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
  MultiscaleReactionProcess() {}
  virtual ~MultiscaleReactionProcess() {}
  virtual void initializeThird()
    { 
      //Set up the following:
      //theMultiscale = isMultiscale species 
      //M = a species on multiscale (isOnMultiscale)
      //N = a normal species that can bind with theMultiscaleSpecies to become
      //    M
      //This must be in initializeSecond or later since we need to know
      //if a species is multiscale, which is only set by the
      //CompartmentProcess in initializeFirst:
      if(A->getIsMultiscale())
        {
          theMultiscale = A;
          N = B;
        }
      else if(B->getIsMultiscale())
        {
          theMultiscale = B;
          N = A;
        }
      else
        {
          THROW_EXCEPTION(ValueError, String(
             getPropertyInterface().getClassName()) + " [" + 
              getFullID().asString() + "]: This process must have at least " +
             "one multiscale substrate species.");
        }
      if(N->getVacantSpecies() == theMultiscale)
        {
          THROW_EXCEPTION(ValueError, String(
             getPropertyInterface().getClassName()) + " [" + 
              getFullID().asString() + "]: The substrate " + 
              getIDString(N) + "'s vacant species is " +
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
      if(C->getIsMultiscale() && D->getIsOnMultiscale())
        {
          M = D;
        }
      else if(C->getIsOnMultiscale() && D->getIsMultiscale())
        {
          M = C;
        }
      else
        {
          THROW_EXCEPTION(ValueError, String(
             getPropertyInterface().getClassName()) + " [" + 
              getFullID().asString() + "]: This process must have at least " +
             "one product species on multiscale.");
        }
      //This must be set in
      //initializeThird since it requires vacant species properties
      //set by DiffusionProcess in initializeSecond:

      //If it is a dissociation reaction,
      //M diffuses on theMultiscale,
      //M unbinds from theMultiscale to become N:
      theMultiscale->setMultiscaleBindIDs(N->getID(), M->getID());
      theMultiscale->setMultiscaleUnbindIDs(M->getID(), N->getID());
      theMultiscale->setDiffusionInfluencedReaction(
            dynamic_cast<DiffusionInfluencedReactionProcess*>(this),
            N->getID(), 1); 
    }
  virtual void bind(Voxel* aVoxel, const unsigned vacantIdx)
    {
      const unsigned index(aVoxel->idx%theStride);
      M->addMoleculeInMulti(aVoxel, vacantIdx, N->getTag(index));
      N->softRemoveMolecule(index);
    }
  virtual void unbind(Voxel* aVoxel)
    {
      const unsigned index(aVoxel->idx%theStride);
      N->addMoleculeExMulti(aVoxel, M->getTag(index));
      M->softRemoveMolecule(index);
    }
  virtual void finalizeReaction()
    {
      DiffusionInfluencedReactionProcess::finalizeReaction();
    }
protected:
  Species* theMultiscale;
};

#endif /* __MultiscaleReactionProcess_hpp */





