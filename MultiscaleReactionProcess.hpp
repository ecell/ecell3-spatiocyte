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
  virtual void initializeSecond()
    {
      //This must be in initializeSecond since we need to know
      //if a species is multiscale, which is only set by the
      //CompartmentProcess in initializeFirst:
      if(A->getIsMultiscale())
        {
          theMultiscale = A;
          theSubstrate = B;
        }
      else if(B->getIsMultiscale())
        {
          theMultiscale = B;
          theSubstrate = A;
        }
      else
        {
          THROW_EXCEPTION(ValueError, String(
             getPropertyInterface().getClassName()) + " [" + 
              getFullID().asString() + "]: This process must have at least " +
             "one multiscale substrate species.");
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
          theProduct = D;
        }
      else if(!C->getIsMultiscale() && D->getIsMultiscale())
        {
          theProduct = C;
        }
      else
        {
          THROW_EXCEPTION(ValueError, String(
             getPropertyInterface().getClassName()) + " [" + 
              getFullID().asString() + "]: This process must have at least " +
             "one multiscale product species.");
        }
    }
  virtual void initializeThird()
    {
      //This must be set in
      //initializeThird since it requires vacant species properties
      //set by DiffusionProcess in initializeSecond:

      //If it is a dissociation reaction,
      //theSubstrate diffuses on theMultiscale,
      //theSubstrate unbinds from theMultiscale to become theProduct:
      if(theSubstrate->getVacantSpecies() == theMultiscale)
        {
          theMultiscale->setMultiscaleUnbindIDs(theSubstrate->getID(),
                                                theProduct->getID());
        }
      //If it is a association reaction,
      //theProduct diffuses on theMultiscale,
      //theSubstrate binds with theMultiscale to become theProduct:
      else if(theProduct->getVacantSpecies() == theMultiscale)
        {
          theMultiscale->setMultiscaleBindIDs(theSubstrate->getID(),
                                              theProduct->getID());
        }
      theMultiscale->setDiffusionInfluencedReaction(
            dynamic_cast<DiffusionInfluencedReactionProcess*>(this),
            B->getID(), 1); 
    }
  virtual void react(Voxel* aVoxel)
    {
      theProduct->addMolecule(aVoxel);
    }
  virtual void finalizeReaction()
    {
      DiffusionInfluencedReactionProcess::finalizeReaction();
      theSubstrate->updateMoleculeList();
    }
protected:
  Species* theSubstrate;
  Species* theProduct;
  Species* theMultiscale;
};

#endif /* __MultiscaleReactionProcess_hpp */

