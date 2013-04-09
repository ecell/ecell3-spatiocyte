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

#include "MicrotubuleProcess.hpp"

LIBECS_DM_INIT(MicrotubuleProcess, Process); 

unsigned MicrotubuleProcess::getLatticeResizeCoord(unsigned aStartCoord)
{
  const unsigned aSize(CompartmentProcess::getLatticeResizeCoord(aStartCoord));
  theMinusSpecies->resetFixedAdjoins();
  theMinusSpecies->setMoleculeRadius(DiffuseRadius);
  thePlusSpecies->resetFixedAdjoins();
  thePlusSpecies->setMoleculeRadius(DiffuseRadius);
  for(unsigned i(0); i != theKinesinSpecies.size(); ++i)
    {
      theKinesinSpecies[i]->setMoleculeRadius(DiffuseRadius);
    }
  return aSize;
}

void MicrotubuleProcess::setCompartmentDimension()
{
  if(Length)
    {
      Subunits = (unsigned)rint(Length/(2*DiffuseRadius));
    }
  Length = Subunits*2*DiffuseRadius;
  Width = Radius*2;
  Height = Radius*2;
  theDimension = 1;
  nLength = Length/(VoxelRadius*2);
  nWidth = Width/(VoxelRadius*2);
  nHeight = Height/(VoxelRadius*2);
  theComp->lengthX = nHeight;
  theComp->lengthY = nWidth;
  theComp->lengthZ = nLength;
}


void MicrotubuleProcess::initializeThird()
{
  if(!isCompartmentalized)
    {
      thePoints.resize(endCoord-subStartCoord);
      initializeVectors();
      initializeFilaments(subunitStart, Filaments, Subunits, nMonomerPitch,
                          theMinusSpecies, subStartCoord);
      elongateFilaments(theVacantSpecies, subStartCoord, Filaments, Subunits,
                        nDiffuseRadius);
      connectFilaments(subStartCoord, Filaments, Subunits);
      interfaceSubunits();
      isCompartmentalized = true;
    }
  theVacantSpecies->setIsPopulated();
  theInterfaceSpecies->setIsPopulated();
  theMinusSpecies->setIsPopulated();
  thePlusSpecies->setIsPopulated();
}


void MicrotubuleProcess::initializeVectors()
{ 
  //Minus end
  Minus.x = -nLength/2;
  Minus.y = 0;
  Minus.z = 0;
  //Rotated Minus end
  theSpatiocyteStepper->rotateX(RotateX, &Minus, -1);
  theSpatiocyteStepper->rotateY(RotateY, &Minus, -1);
  theSpatiocyteStepper->rotateZ(RotateZ, &Minus, -1);
  add_(Minus, Origin);
  //Direction vector from the Minus end to center
  //Direction vector from the Minus end to center
  lengthVector = sub(Origin, Minus);
  //Make direction vector a unit vector
  norm_(lengthVector);
  //Rotated Plus end
  Plus = disp(Minus, lengthVector, nLength);
  setSubunitStart();
}

void MicrotubuleProcess::setSubunitStart()
{
  Point R; //Initialize a random point on the plane attached at the minus end
  if(Minus.x != Plus.x)
    {
      R.y = 10;
      R.z = 30; 
      R.x = (Minus.x*lengthVector.x+Minus.y*lengthVector.y-R.y*lengthVector.y+
             Minus.z*lengthVector.z-R.z*lengthVector.z)/lengthVector.x;
    }
  else if(Minus.y != Plus.y)
    {
      R.x = 10; 
      R.z = 30;
      R.y = (Minus.x*lengthVector.x-R.x*lengthVector.x+Minus.y*lengthVector.y+
             Minus.z*lengthVector.z-R.z*lengthVector.z)/lengthVector.y;
    }
  else
    {
      R.x = 10; 
      R.y = 30;
      R.z = (Minus.x*lengthVector.x-R.x*lengthVector.x+Minus.y*lengthVector.y-
             R.y*lengthVector.y+Minus.z*lengthVector.z)/lengthVector.z;
    }
  //The direction vector from the minus end to the random point, R
  Point D(sub(R, Minus));
  norm_(D);
  subunitStart = disp(Minus, D, nRadius);
}

void MicrotubuleProcess::initializeFilaments(Point& aStartPoint, unsigned aRows,
                                             unsigned aCols, double aRadius,
                                             Species* aVacant,
                                             unsigned aStartCoord)
{
  addCompVoxel(0, 0, aStartPoint, aVacant, aStartCoord, aCols);
  Point U(aStartPoint);
  for(unsigned i(1); i < aRows; ++i)
    {
      double angle(2*M_PI/aRows);
      rotatePointAlongVector(U, Minus, lengthVector, angle);
      disp_(U, lengthVector, aRadius/(aRows-1));
      addCompVoxel(i, 0, U, aVacant, aStartCoord, aCols);
    }
}

// y:width:rows:filaments
// z:length:cols:subunits
void MicrotubuleProcess::connectFilaments(unsigned aStartCoord,
                                          unsigned aRows, unsigned aCols)
{
  for(unsigned i(0); i != aCols; ++i)
    {
      for(unsigned j(0); j != aRows; ++j)
        { 
          if(i > 0)
            { 
              //NORTH-SOUTH
              unsigned a(aStartCoord+j*aCols+i);
              unsigned b(aStartCoord+j*aCols+(i-1));
              connectSubunit(a, b, NORTH, SOUTH);
            }
          else if(Periodic)
            {
              //periodic NORTH-SOUTH
              unsigned a(aStartCoord+j*aCols); 
              unsigned b(aStartCoord+j*aCols+aCols-1);
              connectSubunit(a, b, NORTH, SOUTH);
            }
        }
    }
}

void MicrotubuleProcess::elongateFilaments(Species* aVacant,
                                           unsigned aStartCoord,
                                           unsigned aRows,
                                           unsigned aCols,
                                           double aRadius)
{
  for(unsigned i(0); i != aRows; ++i)
    {
      Voxel* startVoxel(&(*theLattice)[aStartCoord+i*aCols]);
      Point A(*startVoxel->point);
      for(unsigned j(1); j != aCols; ++j)
        {
          disp_(A, lengthVector, aRadius*2);
          if(j != aCols-1)
            {
              addCompVoxel(i, j, A, aVacant, aStartCoord, aCols);
            }
          else
            {
              addCompVoxel(i, j, A, thePlusSpecies, aStartCoord, aCols);
            }
        }
    }
}

//Is inside the parent compartment and confined by the length of the MT:
bool MicrotubuleProcess::isInside(Point& aPoint)
{
  double dist(point2planeDist(aPoint, lengthVector, dot(lengthVector, Minus)));
  if(dist > nDiffuseRadius)
    { 
      dist = point2planeDist(aPoint, lengthVector, dot(lengthVector, Plus));
      if(dist < -nDiffuseRadius)
        {
          return true;
        }
    }
  return false;
}

void MicrotubuleProcess::addNonIntersectInterfaceVoxel(Voxel& aVoxel,
                                                       Point& aPoint)
{
  //Get the distance from the voxel to the center line of the MT:
  double distA(point2lineDist(aPoint, lengthVector, Minus));
  for(unsigned i(0); i != theAdjoiningCoordSize; ++i)
    {
      Voxel& adjoin((*theLattice)[aVoxel.adjoiningCoords[i]]);
      //if(getID(adjoin) != theInterfaceSpecies->getID())
      if(theSpecies[getID(adjoin)]->getIsCompVacant())
        {
          Point pointB(theSpatiocyteStepper->coord2point(adjoin.coord));
          //Get the distance from the adjoin to the center line of the MT:
          double distB(point2lineDist(pointB, lengthVector, Minus));
          //If not on the same side of the MT surface:
          if((distA < nRadius) != (distB < nRadius))
            {
              //If the voxel is nearer to the MT surface:
              if(abs(distA-nRadius) < abs(distB-nRadius))
                {
                  theInterfaceSpecies->addMolecule(&aVoxel);
                  return;
                }
              //If the adjoin is nearer to the MT surface:
              else
                {
                  theInterfaceSpecies->addMolecule(&adjoin);
                }
            }
        }
    }
}

