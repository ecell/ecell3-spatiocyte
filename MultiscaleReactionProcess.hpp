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
  typedef void (MultiscaleReactionProcess::*Method)(Voxel*, Voxel*,
                                              const unsigned, const unsigned);
public:
  LIBECS_DM_OBJECT(MultiscaleReactionProcess, Process)
    {
      INHERIT_PROPERTIES(DiffusionInfluencedReactionProcess);
    }
  MultiscaleReactionProcess():
    isReactInMultiscaleComp(false),
    isReactWithMultiscaleComp(false),
    M(NULL),
    N(NULL),
    M_p(NULL),
    N_p(NULL) {}
  virtual ~MultiscaleReactionProcess() {}
  virtual void initializeThird();
  virtual void initializeMultiscaleWalkBindUnbind();
  virtual void initializeMultiscaleCompReaction();
  virtual void setReactMethod();
  virtual void reactWithMultiscaleComp(Voxel*, Voxel*, const unsigned,
                                       const unsigned);
  virtual void reactInMultiscaleComp(Voxel*, Voxel*, const unsigned,
                                       const unsigned);
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
  virtual void react(Voxel* molA, Voxel* molB, const unsigned indexA,
                     const unsigned indexB)
    {
      (this->*reactM)(molA, molB, indexA, indexB);
    }
protected:
  unsigned getIdx(Species*, Voxel*, const unsigned);
  void reactMuAtoMuC(Voxel*, Voxel*, const unsigned, const unsigned);
  void reactMuBtoMuC(Voxel*, Voxel*, const unsigned, const unsigned);
  void reactAtoC_MuBtoMuD(Voxel*, Voxel*, const unsigned, const unsigned);
  void reactMuAtoMuC_BtoD(Voxel*, Voxel*, const unsigned, const unsigned);
  void reactBtoC_MuAtoMuD(Voxel*, Voxel*, const unsigned, const unsigned);
  void reactMuBtoMuC_AtoD(Voxel*, Voxel*, const unsigned, const unsigned);
  void reactAeqC_MuBtoMuD(Voxel*, Voxel*, const unsigned, const unsigned);
  void reactMuAeqMuC_BtoD(Voxel*, Voxel*, const unsigned, const unsigned);
  void reactBeqC_MuAtoMuD(Voxel*, Voxel*, const unsigned, const unsigned);
  void reactMuBeqMuC_AtoD(Voxel*, Voxel*, const unsigned, const unsigned);
  void reactMuBtoMuC_AeqD(Voxel*, Voxel*, const unsigned, const unsigned);
  void reactBtoC_MuAeqMuD(Voxel*, Voxel*, const unsigned, const unsigned);
  void reactMuAtoMuC_BeqD(Voxel*, Voxel*, const unsigned, const unsigned);
  void reactAtoC_MuBeqMuD(Voxel*, Voxel*, const unsigned, const unsigned);
  void reactAtoC_Multi(Voxel*, Voxel*, const unsigned, const unsigned);
  void reactBtoC_Multi(Voxel*, Voxel*, const unsigned, const unsigned);
  void reactMuAtoMuC_MuBtoMuD(Voxel*, Voxel*, const unsigned, const unsigned);
  void reactMuBeqMuC_MuAtoMuD(Voxel*, Voxel*, const unsigned, const unsigned);
  void throwException(String);
  virtual bool getIsReactWithMultiscaleComp()
    {
      return isReactWithMultiscaleComp;
    }
  virtual bool getIsReactInMultiscaleComp()
    {
      return isReactInMultiscaleComp;
    }
protected:
  bool isReactInMultiscaleComp;
  bool isReactWithMultiscaleComp;
  Species* M;
  Species* N;
  Species* M_p;
  Species* N_p;
  Species* theMultiscale;
  Method reactM;
};

#endif /* __MultiscaleReactionProcess_hpp */





