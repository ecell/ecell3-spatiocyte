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
  Species* createSpecies(System*, String);
  Variable* createVariable(System*, String);
  Species* getSpecies(Variable*);
  std::vector<Species*> getSpecies();
  Point coord2point(unsigned);
  void optimizeSurfaceVoxel(unsigned, Comp*);
  void setSurfaceSubunit(unsigned, Comp*);
  Species* id2species(unsigned short);
  Comp* id2Comp(unsigned short);
  void coord2global(unsigned, unsigned&, unsigned&, unsigned&);
  void point2global(Point, unsigned&, unsigned&, unsigned&);
  Comp* system2Comp(System*);
  bool isBoundaryMol(unsigned, unsigned);
  unsigned getPeriodicMol(unsigned, unsigned, Origin*);
  unsigned global2coord(unsigned, unsigned, unsigned);
  Point getPeriodicPoint(unsigned, unsigned, Origin*);
  void checkLattice();
  void setPeriodicEdge();
  void reset(int);
  unsigned getRowSize();
  unsigned getLayerSize();
  unsigned getColSize();
  unsigned getLatticeSize();
  unsigned short getNullID();
  Point getCenterPoint();
  double getNormalizedVoxelRadius();
  unsigned point2coord(Point&);
  std::vector<Comp*> const& getComps() const;
  Species* variable2species(Variable*);
  void rotateX(double, Point*, int sign=1);
  void rotateY(double, Point*, int sign=1);
  void rotateZ(double, Point*, int sign=1);
  bool isPeriodicEdgeMol(unsigned, Comp*);
  bool isRemovableEdgeMol(unsigned, Comp*);
  double getRowLength();
  double getColLength();
  double getLayerLength();
  double getMinLatticeSpace();
  void updateSpecies();
  void finalizeSpecies();
  unsigned getStartMol();

  virtual GET_METHOD(Real, TimeScale)
  {
      return 0.0;
  }

private:
  void setCompsCenterPoint();
  void setIntersectingCompartmentList();
  void setIntersectingParent();
  void setIntersectingPeers();
  void printProcessParameters();
  void checkSurfaceComp();
  void shuffleAdjoins();
  void setLatticeProperties();
  void checkModel();
  void resizeProcessLattice();
  void initPriorityQueue();
  void initializeFirst();
  void initializeSecond();
  void initializeThird();
  void initializeFourth();
  void initializeFifth();
  void initializeLastOnce();
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
  void populateSpeciesDense(std::vector<Species*>&, unsigned, unsigned);
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
  void concatenateVoxel(Voxel&, unsigned, unsigned, unsigned);
  void concatenateLayers(Voxel&, unsigned, unsigned, unsigned, unsigned);
  void concatenateRows(Voxel&, unsigned, unsigned, unsigned, unsigned);
  void concatenateCols(Voxel&, unsigned, unsigned, unsigned, unsigned);
  void replaceVoxel(unsigned, unsigned);
  void replaceUniVoxel(unsigned, unsigned);
  void setMinMaxSurfaceDimensions(unsigned, Comp*);
  bool isInsideMol(unsigned, Comp*, double);
  bool isSurfaceVoxel(Voxel&, unsigned, Comp*);
  bool isLineVoxel(Voxel&, unsigned, Comp*);
  bool isEnclosedSurfaceVoxel(Voxel&, unsigned, Comp*);
  bool isEnclosedRootSurfaceVoxel(Voxel&, unsigned, Comp*, Comp*);
  bool isPeerMol(unsigned, Comp*);
  bool isLowerPeerMol(unsigned, Comp*);
  bool isRootSurfaceVoxel(Voxel&, unsigned, Comp*);
  bool isParentSurfaceVoxel(Voxel&, unsigned, Comp*);
  bool compartmentalizeVoxel(unsigned, Comp*);
  double getCuboidSpecArea(Comp*);
  unsigned coord2row(unsigned);
  unsigned coord2col(unsigned);
  unsigned coord2layer(unsigned);
  Comp* registerComp(System*, std::vector<Comp*>*);
  Variable* getVariable(System*, String const&);
private:
  bool isInitialized;
  bool isPeriodicEdge;
  bool SearchVacant;
  unsigned short theNullID;
  unsigned LatticeType; 
  unsigned theAdjoinSize;
  unsigned theCellShape;
  unsigned theRowSize;
  unsigned theColSize;
  unsigned theLayerSize;
  unsigned theBioSpeciesSize;
  double VoxelRadius; //r_v
  double theNormalizedVoxelRadius;
  double theHCPl;
  double theHCPx;
  double theHCPy;
  unsigned theNullMol;
  Point theCenterPoint;
  ProcessPriorityQueue thePriorityQueue; 
  std::vector<Species*>::iterator variable2ispecies(Variable*);
  std::vector<Species*> theSpecies;
  std::vector<Comp*> theComps;
  std::vector<Voxel> theLattice;
  std::vector<Process*> theExternInterruptedProcesses;
  std::vector<std::vector<unsigned> > theCoordMols;
};

#endif /* __SpatiocyteStepper_hpp */

