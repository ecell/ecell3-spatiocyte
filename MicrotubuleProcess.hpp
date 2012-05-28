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
      PROPERTYSLOT_SET_GET(Real, OriginX);
      PROPERTYSLOT_SET_GET(Real, OriginY);
      PROPERTYSLOT_SET_GET(Real, OriginZ);
      PROPERTYSLOT_SET_GET(Real, RotateX);
      PROPERTYSLOT_SET_GET(Real, RotateY);
      PROPERTYSLOT_SET_GET(Real, RotateZ);
    }
  MicrotubuleProcess():
    DimerPitch(8e-9),
    Length(100e-9),
    MonomerPitch(4e-9),
    Protofilaments(13),
    Radius(12.5e-9),
    OriginX(0),
    OriginY(0),
    OriginZ(0),
    RotateX(0),
    RotateY(0),
    RotateZ(0),
    theVacantSpecies(NULL) {}
  virtual ~MicrotubuleProcess() {}
  SIMPLE_SET_GET_METHOD(Real, DimerPitch);
  SIMPLE_SET_GET_METHOD(Real, Length);
  SIMPLE_SET_GET_METHOD(Real, MonomerPitch);
  SIMPLE_SET_GET_METHOD(Real, Protofilaments);
  SIMPLE_SET_GET_METHOD(Real, Radius);
  SIMPLE_SET_GET_METHOD(Real, OriginX);
  SIMPLE_SET_GET_METHOD(Real, OriginY);
  SIMPLE_SET_GET_METHOD(Real, OriginZ);
  SIMPLE_SET_GET_METHOD(Real, RotateX);
  SIMPLE_SET_GET_METHOD(Real, RotateY);
  SIMPLE_SET_GET_METHOD(Real, RotateZ);
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
              if(theVacantSpecies)
                {
                  THROW_EXCEPTION(ValueError, String(
                                  getPropertyInterface().getClassName()) +
                                  "[" + getFullID().asString() + 
                                  "]: A MicrotubuleProcess requires only one " +
                                  "vacant variable reference with negative " +
                                  "coefficient as the vacant species of the " +
                                  "microtubule compartment, but " +
                                  getIDString(theVacantSpecies) + " and " +
                                  getIDString(aSpecies) + " are given."); 
                }
              theVacantSpecies = aSpecies;
              theVacantSpecies->setIsVacant();
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
    }

  virtual void initializeThird();
  void initProtofilaments();
  void elongateProtofilaments();
  void connectProtofilaments();
  double getWestPlaneDist(Voxel*);
  void initializeDirectionVector();
  Voxel* getNeighbor(Voxel*, Point&, Voxel*, double&);
  bool notNeighbor(Voxel*, Voxel*);
  bool notShared(Voxel*, Point, Voxel*);
  bool isInsidePlane(Voxel*, Point&);
  bool isLine(Voxel*, Point&);
  bool checkStartVoxel(Voxel*);
  void addVacantVoxel(unsigned int, unsigned int, Point&);
  void removeVacantVoxels(unsigned int);
  void rotatePointAlongVector(Point&, double);
  void connectLatticeVoxels();
  void enlistLatticeVoxels();
  bool addLatticeVoxel(Voxel*, Voxel*);
  bool isValidVoxel(Voxel*);
  bool isInsideVoxel(Voxel*);
  void removeInsideAdjoins(Voxel*);
  void connectLatticeVoxel(Voxel*, Voxel*);
  void connectNorthSouth(unsigned int, unsigned int);
  void connectEastWest(unsigned int, unsigned int);
  void connectSeamEastWest(unsigned int);
  void connectNwSw(unsigned int);
  bool notInNeighbors(Voxel*, Point&);
protected:
  double DimerPitch;
  double Length;
  double MonomerPitch;
  double Protofilaments;
  double Radius;
  double VoxelDiameter;
  double OriginX;
  double OriginY;
  double OriginZ;
  double RotateX;
  double RotateY;
  double RotateZ;
  double offLatticeRadius;
  double latticeRadius;
  unsigned int theDimerSize;
  unsigned int theAdjoiningVoxelSize;
  int tempID;
  Comp* theComp;
  Voxel* theNullVoxel;
  Point T; //Direction vector along the MT axis from Minus to Plus end
  Point M; //Minus end
  Point P; //Plus end
  Point C; //Center point
  std::vector<Voxel*> startVoxels;
  std::vector<std::vector<Voxel*> > vacantVoxels;
  Species* theVacantSpecies;
  std::vector<Voxel> theLattice;
  std::vector<Point> thePoints;
  std::vector<Species*> theKinesinSpecies;
  std::vector<Voxel*> latticeVoxels;
};

#endif /* __MicrotubuleProcess_hpp */




