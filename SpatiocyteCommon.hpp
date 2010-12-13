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


#ifndef __SpatiocyteCommon_hpp
#define __SpatiocyteCommon_hpp

#include <iostream>
#include <iomanip>
#include <math.h>
#include <libecs.hpp>
#include <FullID.hpp>
#include "PriorityQueue.hpp"

USE_LIBECS;
using namespace std;
class SpatiocyteProcess;
class Species;
struct Subunit;
typedef PriorityQueue<SpatiocyteProcess*> ProcessPriorityQueue;
typedef ProcessPriorityQueue::ID ProcessID;

//Compartment type:
#define VOLUME  0
#define SURFACE 1

//Compartment shape:
#define SPHERICAL     0
#define ROD           1
#define CUBIC         2
#define CUBOID        3
#define ELLIPSOID     4
#define CYLINDRICAL   5

//CUBIC shaped compartment surface boundary conditions:
#define REFLECTIVE 0 
#define PERIODIC   1
#define UNIPERIODIC   2

//Number of neighbor voxels for a voxel in the hexagonal close-packed (HCP) lattice: 
#define ADJOINING_VOXEL_SIZE 12

//The 12 adjoining voxels of a voxel in the HCP lattice:
#define NORTH    0 
#define SOUTH    1
#define EAST     2 
#define WEST     3 
#define NW       4 
#define SW       5
#define NE       6
#define SE       7
#define DORSALN  8
#define DORSALS  9
#define VENTRALN 10 
#define VENTRALS 11

#define INNER     0
#define OUTER     1
#define IMMEDIATE 2
#define EXTENDED  3
#define SHARED    4

//Polymerization parameters
#define LARGE_DISTANCE 50
#define MAX_MONOMER_OVERLAP 0.2
#define MAX_IMMEDIATE_DISTANCE 0.2 
#define BIG_NUMBER 1e+20

struct Voxel
{
  //We use short here to maintain the size of Voxel as 128 bytes which helps
  //prefetching. Species ID:
  unsigned short id;
  //Try to limit the adjoiningSize <= 6:
  unsigned short adjoiningSize;
  unsigned int coord;
  Voxel* adjoiningVoxels[ADJOINING_VOXEL_SIZE];
  Subunit* subunit;
  //Contains adjoining and extended surface voxels:
  vector<vector<Voxel*> >* surfaceVoxels;
};

struct Point 
{
  double x;
  double y;
  double z;
};

struct Compartment
{
  bool isEnclosed;
  bool isSurface;
  unsigned short vacantID; 
  int shape;
  double lengthX;
  double lengthY;
  double lengthZ;
  double specVolume;
  double specArea;
  double actualVolume;
  double actualArea;
  System* system;
  Compartment* surfaceSub;
  Point centerPoint;
  Point eastPoint;
  Point westPoint;
  vector<Compartment*> allSubs;
  vector<Compartment*> immediateSubs;
  vector<Compartment*> reactiveComps;
  vector<Species*> species;
  vector<unsigned int> coords;
};

struct Origin
{
  Point point;
  int row;
  int layer;
  int col;
};

struct Bend
{
  double angle;
  double dcm[9];
  double sphereDcm[9];
  double cylinderDcm[9];
};

struct Subunit
{
  unsigned int bendSize;
  //point is the surface point of the voxel:
  Point surfacePoint;
  //subunitPoint is the continuous point of the subunit. 
  //subunitPoint = surfacePoint if the voxel is the origin of the polymer:
  Point subunitPoint;
  //voxel is the actual voxel occupied by the molecule. A shared lipid voxel's 
  //subunit also points to the voxel occupied by the molecule. We need to 
  //differentiate the shared and actual voxels:
  Voxel* voxel;
  Species* species;
  vector<Voxel*> targetVoxels;
  vector<bool> boundBends;
  vector<Voxel*> sourceVoxels;
  vector<Voxel*> sharedLipids;
  vector<Voxel*> tmpVoxels;
  vector<Point> targetPoints;
  vector<Bend> targetBends;
  //contPoints are all the continuous points that are represented by the voxel
  //pointed by this subunit:
  //There shouldn't be duplicate contPoints
  vector<Point> contPoints;
  //the contPointSize is the number of times the same continuous point is used.
  //so there can be duplicates of contPoints
  vector<int> contPointSize;
};

#endif /* __SpatiocyteCommon_hpp */
