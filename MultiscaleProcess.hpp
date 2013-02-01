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


#ifndef __MultiscaleProcess_hpp
#define __MultiscaleProcess_hpp

#include <sstream>
#include "SpatiocyteProcess.hpp"
#include "SpatiocyteSpecies.hpp"

LIBECS_DM_CLASS(MultiscaleProcess, SpatiocyteProcess)
{ 
public:
  LIBECS_DM_OBJECT(MultiscaleProcess, Process)
    {
      INHERIT_PROPERTIES(Process);
    }
  MultiscaleProcess() {}
  virtual ~MultiscaleProcess() {}
  virtual void initialize()
    {
      if(isInitialized)
        {
          return;
        }
      SpatiocyteProcess::initialize();
      for(VariableReferenceVector::iterator
          i(theVariableReferenceVector.begin());
          i != theVariableReferenceVector.end(); ++i)
        {
          Species* aSpecies(theSpatiocyteStepper->variable2species(
                                   (*i).getVariable())); 
          if((*i).getCoefficient())
            {
              if((*i).getCoefficient() < 0)
                {
                  theSubstrate = aSpecies;
                }
              else
                {
                  theProduct = aSpecies;
                }
            }
          else
            {
              theMultiscale = aSpecies;
            }
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
    }
protected:
  Species* theSubstrate;
  Species* theProduct;
  Species* theMultiscale;
};

#endif /* __MultiscaleProcess_hpp */

