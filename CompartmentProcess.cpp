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
  Comp* aComp(theSpatiocyteStepper->system2Comp(getSuperSystem()));
  aComp->interfaceID = theInterfaceSpecies->getID();
  *theComp = *aComp;
  theVacantSpecies->resetFixedAdjoins();
  theVacantSpecies->setMoleculeRadius(DiffuseRadius);
  if(theLipidSpecies)
    {
      theLipidSpecies->resetFixedAdjoins();
      theLipidSpecies->setMoleculeRadius(LipidRadius);
    }
  //The compartment center point (origin):
  Origin = aComp->centerPoint;
  Origin.x += OriginX*aComp->lengthX/2;
  Origin.y += OriginY*aComp->lengthY/2;
  Origin.z += OriginZ*aComp->lengthZ/2;
  setCompartmentDimension();
  theComp->dimension = dimension;
  setLipidCompSpeciesProperties();
  setVacantCompSpeciesProperties();
  subStartCoord = aStartCoord;
  lipStartCoord = aStartCoord+Filaments*Subunits;
  endCoord = lipStartCoord+LipidRows*LipidCols;
  return endCoord-aStartCoord;
}

void CompartmentProcess::setVacantCompSpeciesProperties()
{
  for(unsigned i(0); i != theVacantCompSpecies.size(); ++i)
    {
      theVacantCompSpecies[i]->setDimension(dimension);
      theVacantCompSpecies[i]->setMoleculeRadius(SubunitRadius);
      theVacantCompSpecies[i]->setDiffuseRadius(DiffuseRadius);
      if(theLipidSpecies)
        {
          theVacantCompSpecies[i]->setMultiscaleVacantSpecies(theLipidSpecies);
        }
    }
}

int CompartmentProcess::getCoefficient(Species* aSpecies)
{
  for(VariableReferenceVector::iterator i(theVariableReferenceVector.begin());
      i != theVariableReferenceVector.end(); ++i)
    {
      if(aSpecies->getVariable() == (*i).getVariable()) 
        {
          return (*i).getCoefficient();
        }
    }
  return 0;
}

Species* CompartmentProcess::coefficient2species(int aCoeff)
{
  for(VariableReferenceVector::iterator i(theVariableReferenceVector.begin());
      i != theVariableReferenceVector.end(); ++i)
    {
      if((*i).getCoefficient() == aCoeff)
        {
          return theSpatiocyteStepper->variable2species((*i).getVariable());
        }
    }
  return NULL;
}

void CompartmentProcess::setLipidCompSpeciesProperties()
{
  for(unsigned i(0); i != theLipidCompSpecies.size(); ++i)
    {
      theLipidCompSpecies[i]->setDimension(dimension);
      theLipidCompSpecies[i]->setMoleculeRadius(LipidRadius);
    }
}

void CompartmentProcess::updateResizedLattice()
{
  for(unsigned i(0); i != theVacantCompSpecies.size(); ++i)
    {
      //TODO: replace subStartCoord with theVacantSpecies->getCoord(0) 
      //TODO: replace lipStartCoord with theLipidSpecies->getCoord(0) 
      theVacantCompSpecies[i]->setVacStartCoord(subStartCoord, Filaments,
                                                Subunits);
      theVacantCompSpecies[i]->setLipStartCoord(lipStartCoord, LipidRows,
                                                LipidCols);
    }
}

// y:width:rows
// z:length:cols
void CompartmentProcess::setCompartmentDimension()
{
  Point nearest;
  Point farthest;
  getStartVoxelPoint(subunitStart, nearest, farthest);
  double dist(subunitStart.z-nearest.z+nVoxelRadius-nDiffuseRadius);
  subunitStart.z -= int(dist/(nDiffuseRadius*2))*2*nDiffuseRadius;
  dist = subunitStart.y-nearest.y+nVoxelRadius-nDiffuseRadius;
  unsigned cnt(int(dist/(nDiffuseRadius*sqrt(3))));
  subunitStart.y -= cnt*nDiffuseRadius*sqrt(3);
  if(cnt%2 == 1)
    {
      subunitStart.z += nDiffuseRadius;
    }
  if(Autofit)
    {
      Width = (farthest.y-subunitStart.y+nVoxelRadius)*VoxelRadius*2;
      Length = (farthest.z-subunitStart.z+nVoxelRadius)*VoxelRadius*2;
    }
  if(Length)
    {
      Subunits = (unsigned)(Length/(DiffuseRadius*2));
    }
  if(Width)
    {
      Filaments = (unsigned)((Width-2*DiffuseRadius)/
                                 (DiffuseRadius*sqrt(3)))+1;
    }
  if(Periodic && Filaments%2 != 0)
    {
      ++Filaments;
    }
  //Need to use 2.5 here to avoid rounding off error when calculating
  //LipidRows below:
  Width = 2.5*DiffuseRadius+(Filaments-1)*DiffuseRadius*sqrt(3); 
  if(theLipidSpecies)
    {
      LipidCols = (unsigned)(Length/(LipidRadius*2));
      LipidRows = (unsigned)((Width-2*LipidRadius)/(LipidRadius*sqrt(3)))+1;
    }
  Height = 2*DiffuseRadius;
  if(Filaments == 1)
    {
      dimension = 1;
      Length = Subunits*DiffuseRadius*2;
    }
  else
    {
      dimension = 2;
      //Add DiffuseRadius for the protrusion from hexagonal arrangement:
      Length = Subunits*DiffuseRadius*2+DiffuseRadius;
    }
  //Normalized compartment lengths in terms of lattice voxel radius:
  nLength = Length/(VoxelRadius*2);
  nWidth = Width/(VoxelRadius*2);
  nHeight = Height/(VoxelRadius*2);
  theComp->lengthX = nHeight;
  theComp->lengthY = nWidth;
  theComp->lengthZ = nLength;
  gridCols = (unsigned)rint(nLength/nGridSize);
  gridRows = (unsigned)rint(nWidth/nGridSize);
  theGrid.resize(gridCols*gridRows);
  //Actual surface area = Width*Length
}

void CompartmentProcess::initializeThird()
{
  if(!isCompartmentalized)
    {
      thePoints.resize(endCoord-subStartCoord);
      initializeVectors();
      initializeFilaments(subunitStart, Filaments, Subunits, nDiffuseRadius,
                          theVacantSpecies, subStartCoord);
      elongateFilaments(theVacantSpecies, subStartCoord, Filaments, Subunits,
                        nDiffuseRadius);
      connectFilaments(subStartCoord, Filaments, Subunits);
      interfaceSubunits();
      initializeFilaments(lipidStart, LipidRows, LipidCols, nLipidRadius,
                          theLipidSpecies, lipStartCoord);
      elongateFilaments(theLipidSpecies, lipStartCoord, LipidRows,
                        LipidCols, nLipidRadius);
      connectFilaments(lipStartCoord, LipidRows, LipidCols);
      setDiffuseSize(lipStartCoord, endCoord);
      setSpeciesIntersectLipids();
      isCompartmentalized = true;
    }
  theVacantSpecies->setIsPopulated();
  theLipidSpecies->setIsPopulated();
}

// y:width:rows:filaments
// z:length:cols:subunits
void CompartmentProcess::setSpeciesIntersectLipids()
{
  /*
  for(unsigned i(0); i != theLipidSpecies->size(); ++i)
    {
      Point& aPoint(*(*theLattice)[lipStartCoord+i].point);
      unsigned row((unsigned)((aPoint.y-lipidStart.y)/nGridSize));
      unsigned col((unsigned)((aPoint.z-lipidStart.z)/nGridSize));
      theGrid[col+gridCols*row].push_back(lipStartCoord+i);
    }
    */
  for(unsigned i(0); i != theVacantCompSpecies.size(); ++i)
    {
      theVacantCompSpecies[i]->setIntersectLipids(theLipidSpecies, lipidStart,
                                                  Filaments, Subunits,
                                                  nLipidRadius);
    }
}

void CompartmentProcess::getStartVoxelPoint(Point& start, Point& nearest,
                                            Point& farthest)
{
  Comp* aComp(theSpatiocyteStepper->system2Comp(getSuperSystem()));
  Species* surface(aComp->surfaceSub->vacantSpecies);
  double dist(0);
  Point origin;
  origin.x = 0;
  origin.y = 0;
  origin.z = 0;
  if(surface->size())
    {
      nearest = theSpatiocyteStepper->coord2point(surface->getCoord(0));
      farthest = nearest; 
      origin.x = nearest.x;
      dist = getDistance(&nearest, &origin);
      start = nearest;
    }
  for(unsigned i(1); i < surface->size(); ++i)
    {
      Point aPoint(theSpatiocyteStepper->coord2point(surface->getCoord(i)));
      if(aPoint.x < nearest.x)
        {
          nearest.x = aPoint.x;
          origin.x = aPoint.x;
          dist = getDistance(&aPoint, &origin);
          start = nearest;
        }
      else if(aPoint.x == nearest.x)
        {
          origin.x = aPoint.x;
          double aDist(getDistance(&aPoint, &origin));
          if(aDist < dist)
            {
              dist = aDist;
              start = aPoint;
            }
        }
      if(aPoint.y < nearest.y)
        {
          nearest.y = aPoint.y;
        }
      if(aPoint.z < nearest.z)
        {
          nearest.z = aPoint.z;
        }
      if(aPoint.x > farthest.x)
        {
          farthest.x = aPoint.x;
        }
      if(aPoint.y > farthest.y)
        {
          farthest.y = aPoint.y;
        }
      if(aPoint.z > farthest.z)
        {
          farthest.z = aPoint.z;
        }
    }
}

void CompartmentProcess::initializeVectors()
{
  lengthStart = subunitStart;
  //For Lipid start:
  lengthStart.z -= nDiffuseRadius;
  lengthStart.y -= nDiffuseRadius;

  lengthVector.x = 0;
  lengthVector.y = 0;
  lengthVector.z = 1;
  lengthEnd = disp(lengthStart, lengthVector, nLength);

  widthVector.x = 0;
  widthVector.y = 1;
  widthVector.x = 0;
  widthEnd = disp(lengthEnd, widthVector, nWidth);

  heightVector.x = 1;
  heightVector.y = 0;
  heightVector.z = 0;
  heightEnd = disp(widthEnd, heightVector, nHeight);

  if(theLipidSpecies)
    {
      lipidStart = lengthStart;
      disp_(lipidStart, lengthVector, nLipidRadius);
      disp_(lipidStart, widthVector, nLipidRadius);
    }

  Point center(lengthStart);
  disp_(center, lengthVector, nLength/2);
  disp_(center, widthVector, nWidth/2);
  theComp->centerPoint = center;

  //Set up surface vectors:
  surfaceNormal = cross(lengthVector, widthVector);
  surfaceNormal = norm(surfaceNormal);
  surfaceDisplace = dot(surfaceNormal, widthEnd);
  lengthDisplace = dot(lengthVector, lengthStart);
  lengthDisplaceOpp = dot(lengthVector, lengthEnd);
  widthDisplace = dot(widthVector, widthEnd);
  widthDisplaceOpp = dot(widthVector, lengthEnd);
}

void CompartmentProcess::rotate(Point& V)
{
  theSpatiocyteStepper->rotateX(RotateX, &V, -1);
  theSpatiocyteStepper->rotateY(RotateY, &V, -1);
  theSpatiocyteStepper->rotateZ(RotateZ, &V, -1);
}

void CompartmentProcess::initializeFilaments(Point& aStartPoint, unsigned aRows,
                                             unsigned aCols, double aRadius,
                                             Species* aVacant,
                                             unsigned aStartCoord)
{
  //The first comp voxel must have the aStartCoord:
  addCompVoxel(0, 0, aStartPoint, aVacant, aStartCoord, aCols);
  for(unsigned i(1); i != aRows; ++i)
    {
      Point U(aStartPoint);
      disp_(U, widthVector, i*aRadius*sqrt(3)); 
      if(i%2 == 1)
        {
          disp_(U, lengthVector, -aRadius); 
        }
      addCompVoxel(i, 0, U, aVacant, aStartCoord, aCols);
    }
}

void CompartmentProcess::addCompVoxel(unsigned rowIndex, 
                                      unsigned colIndex,
                                      Point& aPoint,
                                      Species* aVacant,
                                      unsigned aStartCoord,
                                      unsigned aCols)
{
  unsigned aCoord(aStartCoord+rowIndex*aCols+colIndex);
  Voxel& aVoxel((*theLattice)[aCoord]);
  aVoxel.point = &thePoints[aStartCoord-subStartCoord+rowIndex*aCols+colIndex];
  *aVoxel.point = aPoint;
  if(RegularLattice)
    {
      aVoxel.adjoiningSize = theDiffuseSize;
      aVoxel.diffuseSize = theDiffuseSize;
      aVoxel.adjoiningCoords = new unsigned int[theDiffuseSize];
      for(unsigned i(0); i != theDiffuseSize; ++i)
        {
          aVoxel.adjoiningCoords[i] = theNullCoord;
        }
    }
  else
    {
      aVoxel.adjoiningSize = 0;
    }
  aVacant->addCompVoxel(aCoord);
}

void CompartmentProcess::elongateFilaments(Species* aVacant,
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
          addCompVoxel(i, j, A, aVacant, aStartCoord, aCols);
        }
    }
}


void CompartmentProcess::connectSubunit(unsigned a, unsigned b, 
                                        unsigned adjoinA, unsigned adjoinB)
{
  Voxel& voxelA((*theLattice)[a]);
  Voxel& voxelB((*theLattice)[b]);
  if(RegularLattice)
    {
      voxelA.adjoiningCoords[adjoinA] = b;
      voxelB.adjoiningCoords[adjoinB] = a;
    }
  else
    {
      addAdjoin(voxelA, b);
      addAdjoin(voxelB, a);
    }
}

/*
 row0   row1    row2
 fil0   fil1    fil2
 [NW] [ NORTH ] [NE] sub0, col0
      [subunit]      sub1, col1
 [SW] [ SOUTH ] [SE] sub2, col2
 */

// y:width:rows:filaments
// z:length:cols:subunits
void CompartmentProcess::connectFilaments(unsigned aStartCoord,
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
              if(j%2 == 1)
                {
                  if(j+1 < aRows)
                    {
                      //periodic NE-SW 
                      b = aStartCoord+(j+1)*aCols+aCols-1; 
                      connectSubunit(a, b, NE, SW);
                    }
                  else if(j == aRows-1)
                    {
                      //periodic NE-SW 
                      b = aStartCoord+aCols-1; 
                      connectSubunit(a, b, NE, SW);
                    }
                  //periodic NW-SE
                  b = aStartCoord+(j-1)*aCols+aCols-1; 
                  connectSubunit(a, b, NW, SE);
                }
            }
          if(j > 0)
            {
              if(j%2 == 1)
                {
                  //SW-NE
                  unsigned a(aStartCoord+j*aCols+i);
                  unsigned b(aStartCoord+(j-1)*aCols+i); 
                  connectSubunit(a, b, SW, NE);
                  if(i > 0)
                    {
                      //NW-SE
                      b = aStartCoord+(j-1)*aCols+(i-1); 
                      connectSubunit(a, b, NW, SE);
                    }
                }
              else
                {
                  //NW-SE
                  unsigned a(aStartCoord+j*aCols+i);
                  unsigned b(aStartCoord+(j-1)*aCols+i); 
                  connectSubunit(a, b, NW, SE);
                  if(i+1 < aCols)
                    {
                      //SW-NE
                      b = aStartCoord+(j-1)*aCols+(i+1); 
                      connectSubunit(a, b, SW, NE);
                    }
                }
            }
        }
      if(Periodic && aRows > 1)
        { 
          //periodic NW-SE
          unsigned a(aStartCoord+i); //row 0
          unsigned b(aStartCoord+(aRows-1)*aCols+i); 
          connectSubunit(a, b, NW, SE);
          if(i+1 < aCols)
            {
              //periodic SW-NE
              b = aStartCoord+(aRows-1)*aCols+(i+1); 
              connectSubunit(a, b, SW, NE);
            }
        }
    }
}

void CompartmentProcess::setDiffuseSize(unsigned start, unsigned end)
{
  for(unsigned i(start); i != end; ++i)
    {
      Voxel& subunit((*theLattice)[i]);
      subunit.diffuseSize = subunit.adjoiningSize;
    }
}

void CompartmentProcess::addAdjoin(Voxel& aVoxel, unsigned coord)
{
  unsigned* temp(new unsigned[aVoxel.adjoiningSize+1]);
  for(unsigned int i(0); i != aVoxel.adjoiningSize; ++i)
    {
      //Avoid duplicated adjoins:
      if(aVoxel.adjoiningCoords[i] == coord)
        {
          delete[] temp;
          return;
        }
      temp[i] = aVoxel.adjoiningCoords[i];
    }
  delete[] aVoxel.adjoiningCoords;
  temp[aVoxel.adjoiningSize++] = coord;
  aVoxel.adjoiningCoords = temp;
}

void CompartmentProcess::interfaceSubunits()
{
  enlistInterfaceVoxels();
  enlistNonIntersectInterfaceVoxels();
  setDiffuseSize(subStartCoord, lipStartCoord);
  enlistSubunitInterfaceAdjoins();
  theVacantSpecies->setIsPopulated();
  theInterfaceSpecies->setIsPopulated();
}

void CompartmentProcess::enlistInterfaceVoxels()
{
  subunitInterfaces.resize(Filaments*Subunits);
  for(unsigned i(subStartCoord); i != lipStartCoord; ++i)
    {
      Voxel& subunit((*theLattice)[i]);
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
  //Should use SubunitRadius instead of DiffuseRadius since it is the
  //actual size of the subunit:
  if(dist <= (nSubunitRadius+nVoxelRadius)*1.0001) 
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
      //each subunit list has unique (no duplicates) interface voxels:
      subunitInterfaces[subunitCoord-subStartCoord].push_back(voxelCoord);
    }
}

void CompartmentProcess::enlistSubunitInterfaceAdjoins()
{
  for(unsigned i(0); i != subunitInterfaces.size(); ++i)
    {
      for(unsigned j(0); j != subunitInterfaces[i].size(); ++j)
        {
          Voxel& subunit((*theLattice)[i+subStartCoord]);
          Voxel& interface((*theLattice)[subunitInterfaces[i][j]]);
          addAdjoin(interface, i+subStartCoord);
          for(unsigned k(0); k != interface.diffuseSize; ++k)
            {
              unsigned coord(interface.adjoiningCoords[k]);
              Voxel& adjoin((*theLattice)[coord]);
              if(adjoin.id != theInterfaceSpecies->getID())
                {
                  addAdjoin(subunit, coord);
                }
            }
        }
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
  double dist(point2planeDist(aPoint, lengthVector, lengthDisplace));
  if(dist >= 0)
    {
      dist = point2planeDist(aPoint, lengthVector, lengthDisplaceOpp);
      if(dist <= 0)
        {
          dist = point2planeDist(aPoint, widthVector, widthDisplaceOpp);
          if(dist <=0)
            {
              dist = point2planeDist(aPoint, widthVector, widthDisplace);
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

void CompartmentProcess::printParameters()
{
  std::cout << getPropertyInterface().getClassName() << "[" <<
    getFullID().asString() << "]" << std::endl;
  std::cout << "  width:" << Width << " length:" << Length <<
    " area:" << Width*Length << " LipidRows:" << LipidRows << " LipidCols:" <<
    LipidCols << " Filaments:" << Filaments << " Subunits:" << Subunits <<
    std::endl;
  if(theLipidSpecies)
    {
      std::cout << "  " << getIDString(theLipidSpecies) << 
        " number:" << theLipidSpecies->size() << std::endl;
      for(unsigned i(0); i != theLipidCompSpecies.size(); ++i)
        {
          std::cout << "    " << getIDString(theLipidCompSpecies[i]) <<
            " number:" << theLipidCompSpecies[i]->size() << std::endl;
        }
    } 
  std::cout << "  " << getIDString(theVacantSpecies) << 
    " number:" << theVacantSpecies->size() << std::endl;
      for(unsigned i(0); i != theVacantCompSpecies.size(); ++i)
        {
          std::cout << "    " << getIDString(theVacantCompSpecies[i]) <<
            " number:" << theVacantCompSpecies[i]->size() << std::endl;
        }
}
