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
    SearchVacant(false),
    LatticeType(HCP_LATTICE),
    VoxelRadius(10e-9),
    theNormalizedVoxelRadius(0.5) {}
  virtual ~SpatiocyteStepper() {}
  virtual void initialize();
  // need to check interrupt when we suddenly stop the simulation, do we
  // need to update the priority queue?
  virtual void interrupt(Time);
  virtual void step();
  void createSpecies();
  Species* addSpecies(Variable*);
  Species* getSpecies(Variable*);
  std::vector<Species*> getSpecies();
  Point coord2point(unsigned int);
  void optimizeSurfaceVoxel(Voxel&, Comp*);
  void setSurfaceSubunit(Voxel&, Comp*);
  Species* id2species(unsigned short);
  Comp* id2Comp(unsigned short);
  void coord2global(unsigned int, unsigned int&, unsigned int&, unsigned int&);
  void point2global(Point, unsigned int&, unsigned int&, unsigned int&);
  Comp* system2Comp(System*);
  bool isBoundaryCoord(unsigned int, unsigned int);
  unsigned int getPeriodicCoord(unsigned int, unsigned int, Origin*);
  unsigned int global2coord(unsigned int, unsigned int, unsigned int);
  Point getPeriodicPoint(unsigned int, unsigned int, Origin*);
  void checkLattice();
  void setPeriodicEdge();
  void reset(int);
  unsigned int getRowSize();
  unsigned int getLayerSize();
  unsigned int getColSize();
  unsigned int getLatticeSize();
  unsigned short getNullID();
  Point getCenterPoint();
  double getNormalizedVoxelRadius();
  unsigned int point2coord(Point&);
  std::vector<Comp*> const& getComps() const;
  Species* variable2species(Variable*);
  void rotateX(double, Point*, int sign=1);
  void rotateY(double, Point*, int sign=1);
  void rotateZ(double, Point*, int sign=1);
  bool isPeriodicEdgeCoord(unsigned int, Comp*);
  bool isRemovableEdgeCoord(unsigned int, Comp*);
  double getRowLength();
  double getColLength();
  double getLayerLength();
  double getMinLatticeSpace();
  void updateSpecies();
  void finalizeSpecies();
  unsigned int getStartCoord();
private:
  void setCompsCenterPoint();
  void setIntersectingCompartmentList();
  void setIntersectingParent();
  void setIntersectingPeers();
  void printProcessParameters();
  void checkSurfaceComp();
  void shuffleAdjoiningCoords();
  void setLatticeProperties();
  void checkModel();
  void resizeProcessLattice();
  void initPriorityQueue();
  void initProcessSecond();
  void initProcessThird();
  void initProcessFourth();
  void initProcessFifth();
  void initProcessLastOnce();
  void storeSimulationParameters();
  void setSystemSize(System*, double);
  void printSimulationParameters();
  void setCompProperties();
  void initSpecies();
  void readjustSurfaceBoundarySizes();
  void constructLattice();
  void compartmentalizeLattice();
  void concatenatePeriodicSurfaces();
  void registerComps();
  void setCompsProperties();
  void setCompVoxelProperties();
  void populateComps();
  void broadcastLatticeProperties();
  void clearComps();
  void clearComp(Comp*);
  void populateComp(Comp*);
  void populateSpeciesDense(std::vector<Species*>&, unsigned int, unsigned int);
  void populateSpeciesSparse(std::vector<Species*>&);
  void registerCompSpecies(Comp*);
  void setCompProperties(Comp*);
  void setDiffusiveComp(Comp*);
  void setCompCenterPoint(Comp*);
  void setLineVoxelProperties(Comp*);
  void setLineCompProperties(Comp*);
  void setSurfaceVoxelProperties(Comp*);
  void setSurfaceCompProperties(Comp*);
  void setVolumeCompProperties(Comp*);
  void concatenateVoxel(Voxel&, unsigned int, unsigned int, unsigned int);
  void concatenateLayers(Voxel&, unsigned int, unsigned int, unsigned int);
  void concatenateRows(Voxel&, unsigned int, unsigned int, unsigned int);
  void concatenateCols(Voxel&, unsigned int, unsigned int, unsigned int);
  void replaceVoxel(Voxel&, Voxel&);
  void replaceUniVoxel(Voxel&, Voxel&);
  void setMinMaxSurfaceDimensions(unsigned int, Comp*);
  bool isInsideCoord(unsigned int, Comp*, double);
  bool isSurfaceVoxel(Voxel&, Comp*);
  bool isLineVoxel(Voxel&, Comp*);
  bool isEnclosedSurfaceVoxel(Voxel&, Comp*);
  bool isEnclosedRootSurfaceVoxel(Voxel&, Comp*, Comp*);
  bool isPeerCoord(unsigned int, Comp*);
  bool isLowerPeerCoord(unsigned int, Comp*);
  bool isRootSurfaceVoxel(Voxel&, Comp*);
  bool isParentSurfaceVoxel(Voxel&, Comp*);
  bool compartmentalizeVoxel(Voxel&, Comp*);
  double getCuboidSpecArea(Comp*);
  unsigned int coord2row(unsigned int);
  unsigned int coord2col(unsigned int);
  unsigned int coord2layer(unsigned int);
  Comp* registerComp(System*, std::vector<Comp*>*);
  Variable* getVariable(System*, String const&);
private:
  bool isInitialized;
  bool isPeriodicEdge;
  bool SearchVacant;
  unsigned short theNullID;
  unsigned int LatticeType; 
  unsigned int theAdjoiningCoordSize;
  unsigned int theCellShape;
  unsigned int theRowSize;
  unsigned int theColSize;
  unsigned int theLayerSize;
  unsigned int theBioSpeciesSize;
  double VoxelRadius; //r_v
  double theNormalizedVoxelRadius;
  double theHCPl;
  double theHCPx;
  double theHCPy;
  unsigned int theNullCoord;
  Point theCenterPoint;
  ProcessPriorityQueue thePriorityQueue; 
  std::vector<Species*>::iterator variable2ispecies(Variable*);
  std::vector<Species*> theSpecies;
  std::vector<Comp*> theComps;
  std::vector<Voxel> theLattice;
  std::vector<Process*> theQueuedProcesses;
};

#endif /* __SpatiocyteStepper_hpp */

