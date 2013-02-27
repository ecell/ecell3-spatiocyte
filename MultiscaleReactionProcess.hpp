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
      addedVoxels.resize(2);
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
  virtual void bind(Voxel* aVoxel, const unsigned vacantIdx)
    {
      /*
      const unsigned index(aVoxel->idx%theStride);
      if(index < theSubstrates[0]->size())
        {
          removedMols[0].push_back(index);
        }
      aVoxel->idx = addedVoxels[0].size()+theProducts[0]->size()+
        theProducts[0]->getID()*theStride;
      addedVoxels[0].push_back(aVoxel);
      addedIdx.push_back(vacantIdx);
      */
      theSubstrates[0]->softRemoveMolecule(aVoxel);
      theProducts[0]->addMolecule(aVoxel, vacantIdx);
    }
  virtual void unbind(Voxel* aVoxel)
    {
      /*
      const unsigned index(aVoxel->idx%theStride);
      if(index < theSubstrates[1]->size())
        {
          removedMols[1].push_back(index);
        }
      aVoxel->idx = addedVoxels[1].size()+theProducts[1]->size()+
        theProducts[1]->getID()*theStride;
      addedVoxels[1].push_back(aVoxel);
      */
      theSubstrates[1]->softRemoveMolecule(aVoxel);
      theProducts[1]->addMolecule(aVoxel);
    }
  virtual void finalizeReaction()
    {
      /*
      theSubstrates[0]->updateMoleculeList(removedMols[0], addedVoxels[1]);
      removedMols[0].resize(0);
      addedVoxels[1].resize(0);
      theSubstrates[1]->updateMoleculeList(removedMols[1], addedVoxels[0],
                                           addedIdx);
      removedMols[1].resize(0);
      addedVoxels[0].resize(0);
      addedIdx.resize(0);
      */
      DiffusionInfluencedReactionProcess::finalizeReaction();
    }
protected:
  std::vector<Species*> theSubstrates;
  std::vector<Species*> theProducts;
  std::vector<std::vector<Voxel*> > addedVoxels;
  std::vector<std::vector<unsigned> > removedMols;
  std::vector<unsigned> addedIdx;
  Species* theMultiscale;
  double time;
  bool isFinalized;
};

#endif /* __MultiscaleReactionProcess_hpp */





