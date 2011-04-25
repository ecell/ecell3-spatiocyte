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


#ifndef __SpatiocyteStepper_hpp
#define __SpatiocyteStepper_hpp

#include <Stepper.hpp>
#include "SpatiocyteCommon.hpp"

LIBECS_DM_CLASS(SpatiocyteStepper, Stepper)
{ 
public: 
  LIBECS_DM_OBJECT(SpatiocyteStepper, Stepper)
    {
      INHERIT_PROPERTIES(Stepper);
      PROPERTYSLOT_SET_GET(Real, VoxelRadius);
      PROPERTYSLOT_SET_GET(Integer, LatticeType);
      PROPERTYSLOT_SET_GET(Integer, SearchVacant);
    }
  SIMPLE_SET_GET_METHOD(Real, VoxelRadius); 
  SIMPLE_SET_GET_METHOD(Integer, LatticeType); 
  SIMPLE_SET_GET_METHOD(Integer, SearchVacant); 
  SpatiocyteStepper():
    isInitialized(false),
    isPeriodicEdge(false),
    SearchVacant(true),
    LatticeType(HCP_LATTICE),
    VoxelRadius(10e-9),
    theNormalizedVoxelRadius(0.5) {}
  virtual ~SpatiocyteStepper() {}
  virtual void initialize();
  // need to check interrupt when we suddenly stop the simulation, do we
  // need to update the priority queue?
  virtual void interrupt(Time) {}
  virtual void step();
  void createSpecies();
  Species* addSpecies(Variable*);
  Species* getSpecies(Variable*);
  std::vector<Species*> getSpecies();
  Point coord2point(unsigned int);
  void optimizeSurfaceVoxel(Voxel*, Compartment*);
  void setSurfaceSubunit(Voxel*, Compartment*);
  Species* id2species(unsigned short);
  Compartment* id2compartment(unsigned short);
  Voxel* coord2voxel(unsigned int);
  Compartment* system2compartment(System*);
  bool isBoundaryCoord(unsigned int, bool);
  Voxel* getPeriodicVoxel(unsigned int, bool, Origin*);
  Point getPeriodicPoint(unsigned int, bool, Origin*);
  void checkLattice();
  void setPeriodicEdge();
  void reset(int);
  unsigned int getStartCoord();
  unsigned int getRowSize();
  unsigned int getLayerSize();
  unsigned int getColSize();
  unsigned int getLatticeSize();
  Point getCenterPoint();
  double getNormalizedVoxelRadius();
  Voxel* point2voxel(Point);
  std::vector<Compartment*> const& getCompartments() const;
private:
  void setCompartmentsCenterPoint();
  void printProcessParameters();
  void checkSurfaceCompartment();
  void shuffleAdjoiningVoxels();
  void setLatticeProperties();
  void initPriorityQueue();
  void initProcessSecond();
  void initProcessThird();
  void initProcessFourth();
  void initProcessLastOnce();
  void storeSimulationParameters();
  void printSimulationParameters();
  void setSurfaceVoxelProperties();
  void initSpecies();
  void readjustSurfaceBoundarySizes();
  void constructLattice();
  void compartmentalizeLattice();
  void concatenatePeriodicSurfaces();
  void registerCompartments();
  void setCompartmentsProperties();
  void populateCompartments();
  void clearCompartments();
  void clearCompartment(Compartment*);
  void populateCompartment(Compartment*);
  void populateCompartmentUniformly(Compartment*);
  void registerCompartmentSpecies(Compartment*);
  void setCompartmentProperties(Compartment*);
  void removePeriodicEdgeVoxels(Compartment*);
  void removeSurfaces(Compartment*);
  void setReactiveCompartments(Compartment*);
  void setCompartmentCenterPoint(Compartment*);
  void populateCompartmentUniformly(Compartment*, unsigned int);
  void rotateX(double, Point*);
  void rotateY(double, Point*);
  void rotateZ(double, Point*);
  void concatenateVoxel(Voxel*, unsigned int, unsigned int, unsigned int);
  void concatenateLayers(Voxel*, unsigned int, unsigned int, unsigned int);
  void concatenateRows(Voxel*, unsigned int, unsigned int, unsigned int);
  void concatenateCols(Voxel*, unsigned int, unsigned int, unsigned int);
  void coord2global(unsigned int, unsigned int*, unsigned int*, unsigned int*);
  void replaceVoxel(Voxel*, Voxel*);
  void replaceUniVoxel(Voxel*, Voxel*);
  bool isRemovableEdgeCoord(unsigned int, Compartment*);
  bool isReactiveCompartment(Compartment*, Compartment*);
  bool isInsideCoord(unsigned int, Compartment*, double);
  bool isPeriodicEdgeCoord(unsigned int, Compartment*);
  bool isSurfaceVoxel(Voxel*, Compartment*);
  bool isEnclosedSurfaceVoxel(Voxel*, Compartment*);
  bool isPeerVoxel(Voxel*, Compartment*);
  bool compartmentalizeVoxel(Voxel*, Compartment*);
  double getCuboidSpecArea(Compartment*);
  unsigned int coord2row(unsigned int);
  unsigned int coord2col(unsigned int);
  unsigned int coord2layer(unsigned int);
  Compartment* registerCompartment(System*, std::vector<Compartment*>*);
  Variable* getVariable(System*, String const&);
private:
  bool isInitialized;
  bool isPeriodicEdge;
  bool SearchVacant;
  unsigned short theNullID;
  int LatticeType; 
  unsigned int theCellShape;
  unsigned int theStartCoord;
  unsigned int theRowSize;
  unsigned int theColSize;
  unsigned int theLayerSize;
  unsigned int theBioSpeciesSize;
  double VoxelRadius; //r_v
  double theNormalizedVoxelRadius;
  double theHCPk;
  double theHCPh;
  double theHCPl;
  Point theCenterPoint;
  ProcessPriorityQueue thePriorityQueue; 
  std::vector<Species*>::iterator variable2species(Variable*);
  std::vector<Species*> theSpecies;
  std::vector<Compartment*> theCompartments;
  std::vector<Voxel> theLattice;
};

#endif /* __SpatiocyteStepper_hpp */

