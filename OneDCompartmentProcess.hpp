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


#ifndef __OneDCompartmentProcess_hpp
#define __OneDCompartmentProcess_hpp

#include <sstream>
#include <MethodProxy.hpp>
#include "SpatiocyteProcess.hpp"
#include "SpatiocyteSpecies.hpp"

LIBECS_DM_CLASS(OneDCompartmentProcess, SpatiocyteProcess)
{ 
public:
  LIBECS_DM_OBJECT(OneDCompartmentProcess, Process)
    {
      INHERIT_PROPERTIES(Process);
      PROPERTYSLOT_SET_GET(Integer, Quantity);
    }
  OneDCompartmentProcess():
    Quantity(1) {}
  virtual ~OneDCompartmentProcess() {}
  SIMPLE_SET_GET_METHOD(Integer, Quantity);
  virtual void initializeThird();
  virtual void initializeFourth();
  Voxel* getBeginVoxel();
  Voxel* getNeighbor(Voxel*, Point&);
protected:
  Comp* theComp;
  int Quantity;
  Point D; //direction vector from west to east
  Point W; //west point
  Point E; //east point
  std::vector<std::vector<Voxel*> > VacantVoxels;
};

#endif /* __OneDCompartmentProcess_hpp */

