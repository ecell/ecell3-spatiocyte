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


#ifndef __FilamentProcess_hpp
#define __FilamentProcess_hpp

#include <sstream>
#include <MethodProxy.hpp>
#include "SpatiocyteProcess.hpp"
#include "SpatiocyteSpecies.hpp"

LIBECS_DM_CLASS(FilamentProcess, SpatiocyteProcess)
{ 
public:
  LIBECS_DM_OBJECT(FilamentProcess, Process)
    {
      INHERIT_PROPERTIES(Process);
      PROPERTYSLOT_SET_GET(Integer, Filaments);
      PROPERTYSLOT_SET_GET(Integer, Periodic);
      PROPERTYSLOT_SET_GET(Integer, Subunits);
      PROPERTYSLOT_SET_GET(Real, Length);
      PROPERTYSLOT_SET_GET(Real, OriginX);
      PROPERTYSLOT_SET_GET(Real, OriginY);
      PROPERTYSLOT_SET_GET(Real, OriginZ);
      PROPERTYSLOT_SET_GET(Real, RotateX);
      PROPERTYSLOT_SET_GET(Real, RotateY);
      PROPERTYSLOT_SET_GET(Real, RotateZ);
      PROPERTYSLOT_SET_GET(Real, SubunitRadius);
      PROPERTYSLOT_SET_GET(Real, Width);
    }
  FilamentProcess():
    isCompartmentalized(false),
    dimension(1),
    Filaments(1),
    Periodic(0),
    Subunits(1),
    Length(0),
    nVoxelRadius(0.5),
    OriginX(0),
    OriginY(0),
    OriginZ(0),
    RotateX(0),
    RotateY(0),
    RotateZ(0),
    SubunitRadius(0),
    Width(0),
    theVacantSpecies(NULL) {}
  virtual ~FilamentProcess() {}
  SIMPLE_SET_GET_METHOD(Integer, Filaments);
  SIMPLE_SET_GET_METHOD(Integer, Periodic);
  SIMPLE_SET_GET_METHOD(Integer, Subunits);
  SIMPLE_SET_GET_METHOD(Real, Length);
  SIMPLE_SET_GET_METHOD(Real, OriginX);
  SIMPLE_SET_GET_METHOD(Real, OriginY);
  SIMPLE_SET_GET_METHOD(Real, OriginZ);
  SIMPLE_SET_GET_METHOD(Real, RotateX);
  SIMPLE_SET_GET_METHOD(Real, RotateY);
  SIMPLE_SET_GET_METHOD(Real, RotateZ);
  SIMPLE_SET_GET_METHOD(Real, SubunitRadius);
  SIMPLE_SET_GET_METHOD(Real, Width);
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
                                      "]: A FilamentProcess requires only " +
                                      "one vacant variable reference with -1 " +
                                      "coefficient as the vacant species of " +
                                      "the filament compartment, but " +
                                      getIDString(theVacantSpecies) + " and " +
                                      getIDString(aSpecies) + " are given."); 
                    }
                  theVacantSpecies = aSpecies;
                }
            }
          else
            {
              theFilamentSpecies.push_back(aSpecies);
            }
        }
      if(!theFilamentSpecies.size())
        {
          THROW_EXCEPTION(ValueError, String(
                          getPropertyInterface().getClassName()) +
                          "[" + getFullID().asString() + 
                          "]: A FilamentProcess requires at least one " +
                          "nonHD variable reference with zero coefficient " +
                          "as the filament species, but none is given."); 
        }
      if(!theVacantSpecies)
        {
          THROW_EXCEPTION(ValueError, String(
                          getPropertyInterface().getClassName()) +
                          "[" + getFullID().asString() + 
                          "]: A FilamentProcess requires one " +
                          "nonHD variable reference with negative " +
                          "coefficient as the vacant species, " +
                          "but none is given."); 
        }
      if(!SubunitRadius)
        {
          SubunitRadius = theSpatiocyteStepper->getVoxelRadius();
        }
      //Lattice voxel radius:
      VoxelRadius = theSpatiocyteStepper->getVoxelRadius();
      //Normalized off-lattice voxel radius:
      nSubunitRadius = SubunitRadius/(VoxelRadius*2);

      theVacantSpecies->setIsOffLattice();
      for(unsigned i(0); i != theFilamentSpecies.size(); ++i)
        {
          theFilamentSpecies[i]->setIsOffLattice();
        }
    }
  virtual void initializeSecond()
    {
      SpatiocyteProcess::initializeSecond();
      theVacantSpecies->setIsCompVacant();
    }
  virtual unsigned getLatticeResizeCoord(unsigned);
  virtual void initializeThird();
  void addCompVoxel(unsigned, unsigned, Point&);
  void initializeDirectionVector();
  void initializeFilaments();
  void elongateFilaments();
  void connectPeriodic(unsigned);
  void connectNorthSouth(unsigned, unsigned);
  void connectEastWest(unsigned, unsigned);
  void connectSeamEastWest(unsigned);
  void connectNwSw(unsigned);
  void addDirect(Voxel&, unsigned, Voxel&, unsigned);
  void addIndirect(Voxel&, unsigned, Voxel&, unsigned);
  bool initAdjoins(Voxel&);
  void updateAdjoinSize(Voxel&);
  bool inMTCylinder(Point&);
  void rotatePointAlongVector(Point&, double);
  void connectFilaments();
  void addInterfaceVoxel(unsigned, unsigned);
  void setCompartmentDimension();
  void setCompartmentVectors();
  void interfaceSubunits();
  void enlistInterfaceVoxels();
  void enlistNonIntersectInterfaceVoxels();
  void addNonIntersectInterfaceVoxel(Voxel&, Point&);
  bool isInside(Point&);
protected:
  bool isCompartmentalized;
  int tempID;
  unsigned dimension;
  unsigned endCoord;
  unsigned Filaments;
  unsigned Periodic;
  unsigned startCoord;
  unsigned Subunits;
  double filamentDisplace;
  double filamentDisplaceOpp;
  double Length;
  double nLength;
  double nSubunitRadius;
  double nVoxelRadius;
  double nWidth;
  double OriginX;
  double OriginY;
  double OriginZ;
  double RotateX;
  double RotateY;
  double RotateZ;
  double SubunitRadius;
  double subunitDisplace;
  double subunitDisplaceOpp;
  double surfaceDisplace;
  double VoxelRadius;
  double Width;
  Comp* theComp;
  Point T; //Direction vector along the MT axis from Minus to Plus end
  Point M; //Minus end
  Point P; //Plus end
  Point C; //Center point
  Point filamentEnd;
  Point filamentStart;
  Point filamentVector;
  Point subunitVector;
  Point surfaceEnd;
  Point surfaceNormal;
  Species* theVacantSpecies;
  Species* theInterfaceSpecies;
  Variable* theInterfaceVariable;
  std::vector<Point> thePoints;
  std::vector<Species*> theFilamentSpecies;
  std::vector<unsigned> occCoords;
  std::vector<std::vector<unsigned> > subunitInterfaces;
};

#endif /* __FilamentProcess_hpp */




