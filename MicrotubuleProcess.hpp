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


#ifndef __MicrotubuleProcess_hpp
#define __MicrotubuleProcess_hpp

#include <sstream>
#include <MethodProxy.hpp>
#include "SpatiocyteProcess.hpp"
#include "SpatiocyteSpecies.hpp"

LIBECS_DM_CLASS(MicrotubuleProcess, SpatiocyteProcess)
{ 
public:
  LIBECS_DM_OBJECT(MicrotubuleProcess, Process)
    {
      INHERIT_PROPERTIES(Process);
      PROPERTYSLOT_SET_GET(Real, DimerPitch);
      PROPERTYSLOT_SET_GET(Real, Length);
      PROPERTYSLOT_SET_GET(Real, MonomerPitch);
      PROPERTYSLOT_SET_GET(Real, Protofilaments);
      PROPERTYSLOT_SET_GET(Real, Radius);
    }
  MicrotubuleProcess():
    DimerPitch(8e-9),
    Length(100e-9),
    MonomerPitch(4e-9),
    Protofilaments(13),
    Radius(12.5e-9) {}
  virtual ~MicrotubuleProcess() {}
  SIMPLE_SET_GET_METHOD(Real, DimerPitch);
  SIMPLE_SET_GET_METHOD(Real, Length);
  SIMPLE_SET_GET_METHOD(Real, MonomerPitch);
  SIMPLE_SET_GET_METHOD(Real, Protofilaments);
  SIMPLE_SET_GET_METHOD(Real, Radius);
  virtual void initializeThird();
  virtual void initializeFourth();
  void queueStartVoxels();
  double getWestPlaneDist(Voxel*);
  void initializeDirectionVector();
  Voxel* getNeighbor(Voxel*, Point&, Voxel*, double&);
  bool notNeighbor(Voxel*, Voxel*);
  bool notShared(Voxel*, Point, Voxel*);
  bool isInsidePlane(Voxel*, Point&);
  bool isLine(Voxel*, Point&);
  bool checkStartVoxel(Voxel*);
  void addVacantVoxel(unsigned int, Voxel*);
  void removeVacantVoxels(unsigned int);
  void rotatePointAlongVector(Point&, double);
protected:
  double DimerPitch;
  double Length;
  double MonomerPitch;
  double Protofilaments;
  double Radius;
  Comp* theComp;
  Point T; //Direction vector along the MT axis from Minus to Plus end
  Point M; //Minus end
  Point P; //Plus end
  std::vector<Voxel*> startVoxels;
  std::vector<std::vector<Voxel*> > vacantVoxels;
};

#endif /* __MicrotubuleProcess_hpp */


