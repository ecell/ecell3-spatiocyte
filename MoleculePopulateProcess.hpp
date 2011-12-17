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


#ifndef __MoleculePopulateProcess_hpp
#define __MoleculePopulateProcess_hpp

#include "SpatiocyteProcess.hpp"
#include "MoleculePopulateProcessInterface.hpp"

LIBECS_DM_CLASS_EXTRA_1(MoleculePopulateProcess, SpatiocyteProcess, MoleculePopulateProcessInterface)
{ 
public:
  LIBECS_DM_OBJECT(MoleculePopulateProcess, Process)
    {
      INHERIT_PROPERTIES(Process);
      PROPERTYSLOT_SET_GET(Real, OriginX);
      PROPERTYSLOT_SET_GET(Real, OriginY);
      PROPERTYSLOT_SET_GET(Real, OriginZ);
      PROPERTYSLOT_SET_GET(Real, GaussianSigma);
      PROPERTYSLOT_SET_GET(Real, ResetTime);
      PROPERTYSLOT_SET_GET(Real, UniformRadius);
    }
  MoleculePopulateProcess():
    GaussianSigma(0),
    OriginX(0),
    OriginY(0),
    OriginZ(0),
    ResetTime(libecs::INF),
    UniformRadius(1) {}
  virtual ~MoleculePopulateProcess() {}
  SIMPLE_SET_GET_METHOD(Real, OriginX);
  SIMPLE_SET_GET_METHOD(Real, OriginY);
  SIMPLE_SET_GET_METHOD(Real, OriginZ);
  SIMPLE_SET_GET_METHOD(Real, GaussianSigma);
  SIMPLE_SET_GET_METHOD(Real, ResetTime);
  SIMPLE_SET_GET_METHOD(Real, UniformRadius);
  virtual void initializeSecond();
  virtual void populateGaussian(Species*);
  virtual void populateUniformDense(Species*, unsigned int[], unsigned int*);
  virtual void populateUniformSparse(Species* aSpecies);
  virtual void fire();
  virtual void initializeFourth()
    {
      theStepInterval = ResetTime;
      theTime = theStepInterval; 
      thePriorityQueue->move(theQueueID);
    }
protected:
  double GaussianSigma;
  double OriginX;
  double OriginY;
  double OriginZ;
  double ResetTime;
  double UniformRadius;
};

#endif /* __MoleculePopulateProcess_hpp */
