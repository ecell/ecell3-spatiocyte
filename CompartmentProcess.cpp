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

#include "CompartmentProcess.hpp"
#include "Vector.hpp"

LIBECS_DM_INIT(CompartmentProcess, Process); 

unsigned CompartmentProcess::getLatticeResizeCoord(unsigned aStartCoord)
{
  theComp = theSpatiocyteStepper->system2Comp(getSuperSystem());
  theVacantSpecies->resetFixedAdjoins();
  theVacantSpecies->setRadius(SubunitRadius);
  //The compartment center point (origin):
  C = theComp->centerPoint;
  C.x += OriginX*theComp->lengthX/2;
  C.y += OriginY*theComp->lengthY/2;
  C.z += OriginZ*theComp->lengthZ/2;
  setCompartmentDimension();
  for(unsigned i(0); i != theCompartmentSpecies.size(); ++i)
    {
      theCompartmentSpecies[i]->setIsOffLattice();
      theCompartmentSpecies[i]->setDimension(dimension);
      theCompartmentSpecies[i]->setVacantSpecies(theVacantSpecies);
      theCompartmentSpecies[i]->setRadius(SubunitRadius);
    }
  //Normalized compartment lengths in terms of lattice voxel radius:
  nLength = Length/(VoxelRadius*2);
  nWidth = Width/(VoxelRadius*2);
  startCoord = aStartCoord;
  endCoord = startCoord+Filaments*Subunits;
  return endCoord-startCoord;
}

void CompartmentProcess::setCompartmentDimension()
{
  if(Length)
    {
      Subunits = (unsigned)rint(Length/(SubunitRadius*2));
    }
  else
    {
      Length = Subunits*SubunitRadius*2;
    }
  if(Width)
    {
      //Surface with hexagonally packed circles:
      Filaments = (unsigned)rint(Width/(SubunitRadius*sqrt(3)));
      if(Subunits > 1)
        {
          if(Filaments == 1)
            {
              dimension = 1;
            }
          else
            {
              dimension = 2;
            }
        }
      else
        {
          dimension = 0;
        }
    }
  else
    {
      Width = Filaments*SubunitRadius*2;
    }
}

void CompartmentProcess::initializeThird()
{
  if(!isCompartmentalized)
    {
      thePoints.resize(Filaments*Subunits);
      initializeFilaments();
      elongateFilaments();
      connectFilaments();
      setCompartmentVectors();
      interfaceSubunits();
      isCompartmentalized = true;
    }
  theVacantSpecies->setIsPopulated();
}

void CompartmentProcess::initializeFilaments()
{
  initializeDirectionVector();
  Point R; //Initialize a random point on the plane attached at the minus end
  if(M.x != P.x)
    {
      R.y = 10;
      R.z = 30; 
      R.x = (M.x*T.x+M.y*T.y-R.y*T.y+M.z*T.z-R.z*T.z)/T.x;
    }
  else if(M.y != P.y)
    {
      R.x = 10; 
      R.z = 30;
      R.y = (M.x*T.x-R.x*T.x+M.y*T.y+M.z*T.z-R.z*T.z)/T.y;
    }
  else
    {
      R.x = 10; 
      R.y = 30;
      R.z = (M.x*T.x-R.x*T.x+M.y*T.y-R.y*T.y+M.z*T.z)/T.z;
    }
  //The direction vector from the minus end to the random point, R
  Point D(sub(R, M)); 
  D = norm(D);
  //The start point of the first protofilament:
  Point S(M); 
  addCompVoxel(0, 0, S);
  for(unsigned i(1); i != Filaments; ++i)
    {
      Point U(S);
      disp_(U, D, i*nSubunitRadius*sqrt(3)); 
      if(i%2 == 1)
        {
          disp_(U, T, nSubunitRadius); 
        }
      addCompVoxel(i, 0, U);
    }
}

void CompartmentProcess::initializeDirectionVector()
{ 
  /*
   * MEnd = {Mx, My, Mz};(*minus end*) 
   * PEnd = {Px, Py, Pz};(*plus end*)
   * MTAxis = (PEnd - MEnd)/Norm[PEnd - MEnd] (*direction vector along the MT
   * long axis*)
   */
  //Minus end
  M.x = -Length/2;
  M.y = 0;
  M.z = 0;
  //Rotated Minus end
  theSpatiocyteStepper->rotateX(RotateX, &M, -1);
  theSpatiocyteStepper->rotateY(RotateY, &M, -1);
  theSpatiocyteStepper->rotateZ(RotateZ, &M, -1);
  M = add(M, C);
  //Direction vector from the Minus end to center
  T = sub(C, M);
  //Make T a unit vector
  T = norm(T);
  //Rotated Plus end
  P = disp(M, T, Length);
}

void CompartmentProcess::addCompVoxel(unsigned filamentIndex, 
                                   unsigned subunitIndex, Point& aPoint)
{
  unsigned aCoord(startCoord+filamentIndex*Subunits+subunitIndex);
  Voxel& aVoxel((*theLattice)[aCoord]);
  aVoxel.point = &thePoints[filamentIndex*Subunits+subunitIndex];
  *aVoxel.point = aPoint;
  aVoxel.adjoiningCoords = new unsigned[theAdjoiningCoordSize];
  aVoxel.diffuseSize = 2;
  for(unsigned i(0); i != theAdjoiningCoordSize; ++i)
    {
      aVoxel.adjoiningCoords[i] = theNullCoord;
    }
  theVacantSpecies->addCompVoxel(aCoord);
}

void CompartmentProcess::elongateFilaments()
{
  for(unsigned i(0); i != Filaments; ++i)
    {
      Voxel* startVoxel(&(*theLattice)[startCoord+i*Subunits]);
      Point A(*startVoxel->point);
      for(unsigned j(1); j != Subunits; ++j)
        {
          disp_(A, T, nSubunitRadius*2);
          addCompVoxel(i, j, A);
        }
    }
}

void CompartmentProcess::connectFilaments()
{
  for(unsigned i(0); i != Subunits; ++i)
    {
      for(unsigned j(0); j != Filaments; ++j)
        { 
          if(i > 0)
            { 
              connectNorthSouth(i, j);
            }
          else if(Periodic)
            {
              connectPeriodic(j);
            }
          if(j > 0)
            {
              connectEastWest(i, j);
            }
        }
      /*
      if(Filaments > 2)
        {
          connectSeamEastWest(i);
        }
        */
      if(i > 0)
        {
          connectNwSw(i);
        }
    }
}

void CompartmentProcess::connectPeriodic(unsigned j)
{
  unsigned a(startCoord+j*Subunits+Subunits-1);
  Voxel& aVoxel((*theLattice)[a]);
  unsigned b(startCoord+j*Subunits); 
  Voxel& adjoin((*theLattice)[b]);
  aVoxel.adjoiningCoords[NORTH] = b;
  adjoin.adjoiningCoords[SOUTH] = a;
  aVoxel.adjoiningSize = 2;
  adjoin.adjoiningSize = 2;
}

void CompartmentProcess::connectNorthSouth(unsigned i, unsigned j)
{
  unsigned a(startCoord+j*Subunits+(i-1));
  Voxel& aVoxel((*theLattice)[a]);
  unsigned b(startCoord+j*Subunits+i);
  Voxel& adjoin((*theLattice)[b]);
  aVoxel.adjoiningCoords[NORTH] = b;
  adjoin.adjoiningCoords[SOUTH] = a;
  aVoxel.adjoiningSize = 2;
  adjoin.adjoiningSize = 2;
}

void CompartmentProcess::connectEastWest(unsigned i, unsigned j)
{
  unsigned a(startCoord+j*Subunits+i);
  Voxel& aVoxel((*theLattice)[a]);
  unsigned b(startCoord+(j-1)*Subunits+i); 
  Voxel& adjoin((*theLattice)[b]);
  aVoxel.adjoiningCoords[aVoxel.adjoiningSize++] = b;
  adjoin.adjoiningCoords[adjoin.adjoiningSize++] = a;
}

void CompartmentProcess::connectSeamEastWest(unsigned i)
{
  unsigned a(startCoord+i);
  Voxel& aVoxel((*theLattice)[a]);
  unsigned b(startCoord+(Filaments-1)*Subunits+i); 
  Voxel& adjoin((*theLattice)[b]);
  aVoxel.adjoiningCoords[aVoxel.adjoiningSize++] = b;
  adjoin.adjoiningCoords[adjoin.adjoiningSize++] = a;
}

void CompartmentProcess::connectNwSw(unsigned i)
{
  unsigned a(startCoord+i);
  Voxel& aVoxel((*theLattice)[a]);
  unsigned b(startCoord+(Filaments-1)*Subunits+(i-1)); 
  Voxel& adjoin((*theLattice)[b]);
  aVoxel.adjoiningCoords[aVoxel.adjoiningSize++] = b;
  adjoin.adjoiningCoords[adjoin.adjoiningSize++] = a;
}

void CompartmentProcess::setCompartmentVectors()
{
  filamentStart = *(*theLattice)[startCoord].point;
  filamentEnd = *(*theLattice)[startCoord+Subunits-1].point;
  subunitVector = sub(filamentEnd, filamentStart);
  subunitVector = norm(subunitVector);
  filamentStart = disp(filamentStart, subunitVector, -nSubunitRadius);
  filamentEnd = disp(filamentEnd, subunitVector, nSubunitRadius*(1+sqrt(3)));
  if(dimension == 2)
    {
      surfaceEnd = *(*theLattice)[startCoord+
        (Filaments-1)*Subunits+Subunits-1].point;
      surfaceEnd = disp(surfaceEnd, subunitVector, nSubunitRadius*(1+sqrt(3)));
      filamentVector = sub(filamentEnd, surfaceEnd);
      filamentVector = norm(filamentVector);
      surfaceEnd = disp(surfaceEnd, filamentVector, -nSubunitRadius);
      filamentEnd = disp(filamentEnd, filamentVector, nSubunitRadius);
      surfaceNormal = cross(subunitVector, filamentVector);
      surfaceNormal = norm(surfaceNormal);
      surfaceDisplace = dot(surfaceNormal, surfaceEnd);
      subunitDisplace = dot(subunitVector, filamentStart);
      subunitDisplaceOpp = dot(subunitVector, filamentEnd);
      filamentDisplace = dot(filamentVector, surfaceEnd);
      filamentDisplaceOpp = dot(filamentVector, filamentEnd);
    }
}

void CompartmentProcess::interfaceSubunits()
{
  enlistInterfaceVoxels();
  enlistNonIntersectInterfaceVoxels();
  theVacantSpecies->setIsPopulated();
  theInterfaceSpecies->setIsPopulated();
}

void CompartmentProcess::enlistInterfaceVoxels()
{
  subunitInterfaces.resize(Filaments*Subunits);
  for(unsigned i(startCoord); i != endCoord; ++i)
    {
      Voxel& subunit((*theLattice)[i]);
      subunit.diffuseSize = subunit.adjoiningSize;
      Point center(*subunit.point);
      Point bottomLeft(*subunit.point);
      Point topRight(*subunit.point);
      bottomLeft.x -= nSubunitRadius+theSpatiocyteStepper->getColLength();
      bottomLeft.y -= nSubunitRadius+theSpatiocyteStepper->getLayerLength();
      bottomLeft.z -= nSubunitRadius+theSpatiocyteStepper->getRowLength();
      topRight.x += nSubunitRadius+theSpatiocyteStepper->getColLength();
      topRight.y += nSubunitRadius+theSpatiocyteStepper->getLayerLength();
      topRight.z += nSubunitRadius+theSpatiocyteStepper->getRowLength();
      unsigned blRow(0);
      unsigned blLayer(0);
      unsigned blCol(0);
      theSpatiocyteStepper->point2global(bottomLeft, blRow, blLayer, blCol);
      unsigned trRow(0);
      unsigned trLayer(0);
      unsigned trCol(0);
      theSpatiocyteStepper->point2global(topRight, trRow, trLayer, trCol);
      for(unsigned j(blRow); j <= trRow; ++j)
        {
          for(unsigned k(blLayer); k <= trLayer; ++k)
            {
              for(unsigned l(blCol); l <= trCol; ++l)
                {
                  unsigned m(theSpatiocyteStepper->global2coord(j, k, l));
                  addInterfaceVoxel(i, m);
                }
            }
        }
    }
}

void CompartmentProcess::addInterfaceVoxel(unsigned subunitCoord,
                                        unsigned voxelCoord)
{ 
  Voxel& subunit((*theLattice)[subunitCoord]);
  Point subunitPoint(*subunit.point);
  Point voxelPoint(theSpatiocyteStepper->coord2point(voxelCoord));
  double dist(getDistance(&subunitPoint, &voxelPoint));
  if(dist <= nSubunitRadius+nVoxelRadius) 
    {
      Voxel& voxel((*theLattice)[voxelCoord]);
      //theSpecies[6]->addMolecule(&voxel);
      //Insert voxel in the list of interface voxels if was not already:
      //if(voxel.id == theComp->vacantSpecies->getID())
      if(voxel.id != theInterfaceSpecies->getID())
        {
          //theSpecies[5]->addMolecule(&voxel);
          theInterfaceSpecies->addMolecule(&voxel);
        }
      subunitInterfaces[subunitCoord-startCoord].push_back(voxelCoord);
    }
}

void CompartmentProcess::enlistNonIntersectInterfaceVoxels()
{
  for(unsigned i(0); i != theInterfaceSpecies->size(); ++i)
    {
      unsigned voxelCoord(theInterfaceSpecies->getCoord(i));
      Voxel& anInterface((*theLattice)[voxelCoord]);
      for(unsigned j(0); j != theAdjoiningCoordSize; ++j)
        {
          Voxel& adjoin((*theLattice)[anInterface.adjoiningCoords[j]]);
          if(adjoin.id != theInterfaceSpecies->getID())
            {
              Point aPoint(theSpatiocyteStepper->coord2point(adjoin.coord));
              if(isInside(aPoint))
                {
                  addNonIntersectInterfaceVoxel(adjoin, aPoint);
                }
            }
        }
    }
}

bool CompartmentProcess::isInside(Point& aPoint)
{
  double dist(point2planeDist(aPoint, subunitVector, subunitDisplace));
  if(dist >= 0)
    {
      dist = point2planeDist(aPoint, subunitVector, subunitDisplaceOpp);
      if(dist <= 0)
        {
          dist = point2planeDist(aPoint, filamentVector, filamentDisplaceOpp);
          if(dist <=0)
            {
              dist = point2planeDist(aPoint, filamentVector, filamentDisplace);
              if(dist >= 0)
                {
                  return true;
                }
            }
        }
    }
  return false;
}

void CompartmentProcess::addNonIntersectInterfaceVoxel(Voxel& aVoxel,
                                                    Point& aPoint)
{
  double distA(point2planeDist(aPoint, surfaceNormal, surfaceDisplace));
  for(unsigned i(0); i != theAdjoiningCoordSize; ++i)
    {
      Voxel& adjoin((*theLattice)[aVoxel.adjoiningCoords[i]]);
      if(adjoin.id != theInterfaceSpecies->getID())
        {
          Point pointB(theSpatiocyteStepper->coord2point(adjoin.coord));
          double distB(point2planeDist(pointB, surfaceNormal, surfaceDisplace));
          //if not on the same side of the plane:
          if((distA < 0) != (distB < 0))
            {
              if(abs(distA) < abs(distB))
                { 
                  //theSpecies[6]->addMolecule(&aVoxel);
                  theInterfaceSpecies->addMolecule(&aVoxel);
                  return;
                }
              else
                {
                  //theSpecies[6]->addMolecule(&adjoin);
                  theInterfaceSpecies->addMolecule(&adjoin);
                }
            }
        }
    }
}

