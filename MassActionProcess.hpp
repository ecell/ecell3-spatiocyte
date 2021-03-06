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


#ifndef __MassActionProcess_hpp
#define __MassActionProcess_hpp

#include <ContinuousProcess.hpp>
#include "SpatiocyteStepper.hpp"
#include "SpatiocyteCommon.hpp"

LIBECS_DM_CLASS(MassActionProcess, ContinuousProcess)
{ 
public:
  LIBECS_DM_OBJECT(MassActionProcess, Process)
    {
      INHERIT_PROPERTIES(Process);
      PROPERTYSLOT_SET_GET(Real, k);
    }
  MassActionProcess():
    k(0),
    theSpace(0) {}
  virtual ~MassActionProcess() {}
  SIMPLE_SET_GET_METHOD(Real, k);
  virtual void fire();
  virtual void initialize()
    {
      Process::initialize();
      declareUnidirectional();
    }
protected:
  double k;
  double theSpace;
};

#endif /* __MassActionProcess_hpp */

