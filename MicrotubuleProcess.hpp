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
#include "CompartmentProcess.hpp"
#include "SpatiocyteSpecies.hpp"

LIBECS_DM_CLASS(MicrotubuleProcess, CompartmentProcess)
{ 
public:
  LIBECS_DM_OBJECT(MicrotubuleProcess, Process)
    {
      INHERIT_PROPERTIES(CompartmentProcess);
      PROPERTYSLOT_SET_GET(Real, DimerPitch);
      PROPERTYSLOT_SET_GET(Real, Length);
      PROPERTYSLOT_SET_GET(Real, MonomerPitch);
      PROPERTYSLOT_SET_GET(Real, Radius);
    }
  MicrotubuleProcess():
    DimerPitch(8e-9),
    Length(100e-9),
    MonomerPitch(4e-9),
    theMinusSpecies(NULL),
    thePlusSpecies(NULL)
  {
    Filaments = 13;
    Autofit = 0;
    RegularLattice = 0;
  }
  virtual ~MicrotubuleProcess() {}
  SIMPLE_SET_GET_METHOD(Real, DimerPitch);
  SIMPLE_SET_GET_METHOD(Real, Length);
  SIMPLE_SET_GET_METHOD(Real, MonomerPitch);
  SIMPLE_SET_GET_METHOD(Real, Radius);
  virtual void prepreinitialize()
    {
      SpatiocyteProcess::prepreinitialize();
      theInterfaceVariable = createVariable("Interface");
    }
  virtual void initialize()
    {
      if(isInitialized)
        {
          return;
        }
      SpatiocyteProcess::initialize();
      theInterfaceSpecies = theSpatiocyteStepper->addSpecies(
                                                       theInterfaceVariable);
      theInterfaceSpecies->setIsInterface();
      for(VariableReferenceVector::iterator
          i(theVariableReferenceVector.begin());
          i != theVariableReferenceVector.end(); ++i)
        {
          Species* aSpecies(theSpatiocyteStepper->variable2species(
                                   (*i).getVariable())); 
          if((*i).getCoefficient())
            {
              if((*i).getCoefficient() == -1)
                {
                  if(theVacantSpecies)
                    {
                      THROW_EXCEPTION(ValueError, String(
                                      getPropertyInterface().getClassName()) +
                                      "[" + getFullID().asString() + 
                                      "]: A MicrotubuleProcess requires only " +
                                      "one vacant variable reference with -1 " +
                                      "coefficient as the vacant species of " +
                                      "the microtubule compartment, but " +
                                      getIDString(theVacantSpecies) + " and " +
                                      getIDString(aSpecies) + " are given."); 
                    }
                  theVacantSpecies = aSpecies;
                }
              else if((*i).getCoefficient() == -2)
                {
                  if(theMinusSpecies)
                    {
                      THROW_EXCEPTION(ValueError, String(
                                      getPropertyInterface().getClassName()) +
                                      "[" + getFullID().asString() + 
                                      "]: A MicrotubuleProcess requires only " +
                                      "one variable reference with -2 " +
                                      "coefficient as the minus end species " +
                                      "of the microtubule compartment, but " +
                                      getIDString(theMinusSpecies) + " and " +
                                      getIDString(aSpecies) + " are given."); 
                    }
                  theMinusSpecies = aSpecies;
                }
              else if((*i).getCoefficient() == -3)
                {
                  if(thePlusSpecies)
                    {
                      THROW_EXCEPTION(ValueError, String(
                                      getPropertyInterface().getClassName()) +
                                      "[" + getFullID().asString() + 
                                      "]: A MicrotubuleProcess requires only " +
                                      "one variable reference with -3 " +
                                      "coefficient as the plus end species " +
                                      "of the microtubule compartment, but " +
                                      getIDString(thePlusSpecies) + " and " +
                                      getIDString(aSpecies) + " are given."); 
                    }
                  thePlusSpecies = aSpecies;
                }
            }
          else
            {
              theKinesinSpecies.push_back(aSpecies);
            }
        }
      if(!theKinesinSpecies.size())
        {
          THROW_EXCEPTION(ValueError, String(
                          getPropertyInterface().getClassName()) +
                          "[" + getFullID().asString() + 
                          "]: A MicrotubuleProcess requires at least one " +
                          "nonHD variable reference with zero coefficient " +
                          "as the kinesin species, but none is given."); 
        }
      if(!theVacantSpecies)
        {
          THROW_EXCEPTION(ValueError, String(
                          getPropertyInterface().getClassName()) +
                          "[" + getFullID().asString() + 
                          "]: A MicrotubuleProcess requires one " +
                          "nonHD variable reference with negative " +
                          "coefficient as the vacant species, " +
                          "but none is given."); 
        }
      if(!theMinusSpecies)
        {
          theMinusSpecies = theVacantSpecies;
        }
      if(!thePlusSpecies)
        {
          thePlusSpecies = theVacantSpecies;
        }
      if(!DiffuseRadius)
        {
          DiffuseRadius = theSpatiocyteStepper->getVoxelRadius();
        }
      if(!SubunitRadius)
        {
          SubunitRadius = DiffuseRadius;
        }
      VoxelRadius = theSpatiocyteStepper->getVoxelRadius();
      //Normalized off-lattice voxel radius:
      nSubunitRadius = SubunitRadius/(VoxelRadius*2);
      nDiffuseRadius = DiffuseRadius/(VoxelRadius*2);
      nDimerPitch = DimerPitch/(VoxelRadius*2);
      nLength = Length/(VoxelRadius*2);
      nMonomerPitch = MonomerPitch/(VoxelRadius*2);
    }
  virtual void initializeFirst()
    {
      CompartmentProcess::initializeFirst();
      theMinusSpecies->setIsOffLattice();
      theMinusSpecies->setComp(theComp);
      theMinusSpecies->setIsCompVacant();
      thePlusSpecies->setIsOffLattice();
      thePlusSpecies->setComp(theComp);
      thePlusSpecies->setIsCompVacant();
      for(unsigned i(0); i != theKinesinSpecies.size(); ++i)
        {
          theKinesinSpecies[i]->setIsOffLattice();
          theKinesinSpecies[i]->setDimension(1);
          theKinesinSpecies[i]->setVacantSpecies(theVacantSpecies);
          theKinesinSpecies[i]->setComp(theComp);
        }
    }
  virtual unsigned getLatticeResizeCoord(unsigned);
  /*
  virtual void initializeThird();
  void initializeFilaments(Point&, unsigned, unsigned, double, Species*,
                           unsigned);
                           */
protected:
  double DimerPitch;
  double Length;
  double latticeRadius;
  double MonomerPitch;
  double nDimerPitch;
  double nLength;
  double nMonomerPitch;
  double offLatticeRadius;
  double Radius;
  unsigned theDimerSize;
  Point T; //Direction vector along the MT axis from Minus to Plus end
  Point M; //Minus end
  Point P; //Plus end
  Species* theMinusSpecies;
  Species* thePlusSpecies;
  std::vector<Species*> theKinesinSpecies;
  std::vector<unsigned> occCoords;
};

#endif /* __MicrotubuleProcess_hpp */




