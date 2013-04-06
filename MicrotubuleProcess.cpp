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
                          theVacantSpecies, subStartCoord);
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
  add_(M, Origin);
  //Direction vector from the Minus end to center
  lengthVector = sub(Origin, M);
  //Make direction vector a unit vector
  norm_(lengthVector);
  //Rotated Plus end
  P = disp(M, lengthVector, Length);
  setSubunitStart();
}

void MicrotubuleProcess::setSubunitStart()
{
  Point R; //Initialize a random point on the plane attached at the minus end
  if(M.x != P.x)
    {
      R.y = 10;
      R.z = 30; 
      R.x = (M.x*lengthVector.x+M.y*lengthVector.y-R.y*lengthVector.y+
             M.z*lengthVector.z-R.z*lengthVector.z)/lengthVector.x;
    }
  else if(M.y != P.y)
    {
      R.x = 10; 
      R.z = 30;
      R.y = (M.x*lengthVector.x-R.x*lengthVector.x+M.y*lengthVector.y+
             M.z*lengthVector.z-R.z*lengthVector.z)/lengthVector.y;
    }
  else
    {
      R.x = 10; 
      R.y = 30;
      R.z = (M.x*lengthVector.x-R.x*lengthVector.x+M.y*lengthVector.y-
             R.y*lengthVector.y+M.z*lengthVector.z)/lengthVector.z;
    }
  //The direction vector from the minus end to the random point, R
  Point D(sub(R, M));
  norm_(D);
  subunitStart = disp(M, D, nRadius);
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
      rotatePointAlongVector(U, M, lengthVector, angle);
      disp_(U, lengthVector, aRadius/(aRows-1));
      addCompVoxel(i, 0, U, aVacant, aStartCoord, aCols);
    }
}

/*
void MicrotubuleProcess::initializeThird()
{
  if(!isCompartmentalized)
    {
      thePoints.resize(Protofilaments*theDimerSize);
      initializeProtofilaments();
      elongateProtofilaments();
      connectProtofilaments();
      enlistLatticeVoxels();
      isCompartmentalized = true;
    }
  theVacantSpecies->setIsPopulated();
  theMinusSpecies->setIsPopulated();
  thePlusSpecies->setIsPopulated();
}

void MicrotubuleProcess::addCompVoxel(unsigned protoIndex,
                                      unsigned dimerIndex, Point& aPoint)
{
  unsigned aCoord(startCoord+protoIndex*theDimerSize+dimerIndex);
  Voxel& aVoxel((*theLattice)[aCoord]);
  aVoxel.point = &thePoints[protoIndex*theDimerSize+dimerIndex];
  *aVoxel.point = aPoint;
  aVoxel.adjoiningCoords = new unsigned[theAdjoiningCoordSize];
  aVoxel.diffuseSize = 2;
  for(unsigned i(0); i != theAdjoiningCoordSize; ++i)
    {
      aVoxel.adjoiningCoords[i] = theNullCoord;
    }
  if(!dimerIndex)
    {
      theMinusSpecies->addCompVoxel(aCoord);
    }
  else if(dimerIndex == theDimerSize-1)
    { 
      thePlusSpecies->addCompVoxel(aCoord);
    }
  else
    {
      theVacantSpecies->addCompVoxel(aCoord);
    }
}


void MicrotubuleProcess::elongateProtofilaments()
{
  //std::cout << "proto:" << Protofilaments << " dimer:" << theDimerSize << std::endl;
  for(unsigned i(0); i != Protofilaments; ++i)
    {
      Voxel* startVoxel(&(*theLattice)[startCoord+i*theDimerSize]);
      Point A(*startVoxel->point);
      for(unsigned j(1); j != theDimerSize; ++j)
        {
          A.x += DimerPitch*T.x;
          A.y += DimerPitch*T.y;
          A.z += DimerPitch*T.z;
          addCompVoxel(i, j, A);
        }
    }
}

void MicrotubuleProcess::connectProtofilaments()
{
  for(unsigned i(0); i != theDimerSize; ++i)
    {
      for(unsigned j(0); j != Protofilaments; ++j)
        { 
          if(i > 0)
            { 
              connectNorthSouth(i, j);
            }
          else if(Periodic)
            {
              connectPeriodic(j);
            }
          /*
          if(j > 0)
            {
              connectEastWest(i, j);
            }
            *
        }
      /*
      if(Protofilaments > 2)
        {
          connectSeamEastWest(i);
        }
      if(i > 0)
        {
          connectNwSw(i);
        }
        *
    }
}

void MicrotubuleProcess::connectPeriodic(unsigned j)
{
  unsigned a(startCoord+j*theDimerSize+theDimerSize-1);
  Voxel& aVoxel((*theLattice)[a]);
  unsigned b(startCoord+j*theDimerSize); 
  Voxel& adjoin((*theLattice)[b]);
  aVoxel.adjoiningCoords[NORTH] = b;
  adjoin.adjoiningCoords[SOUTH] = a;
  aVoxel.adjoiningSize = 2;
  adjoin.adjoiningSize = 2;
}

void MicrotubuleProcess::connectNorthSouth(unsigned i, unsigned j)
{
  unsigned a(startCoord+j*theDimerSize+(i-1));
  Voxel& aVoxel((*theLattice)[a]);
  unsigned b(startCoord+j*theDimerSize+i);
  Voxel& adjoin((*theLattice)[b]);
  aVoxel.adjoiningCoords[NORTH] = b;
  adjoin.adjoiningCoords[SOUTH] = a;
  aVoxel.adjoiningSize = 2;
  adjoin.adjoiningSize = 2;
}

void MicrotubuleProcess::connectEastWest(unsigned i, unsigned j)
{
  unsigned a(startCoord+j*theDimerSize+i);
  Voxel& aVoxel((*theLattice)[a]);
  unsigned b(startCoord+(j-1)*theDimerSize+i); 
  Voxel& adjoin((*theLattice)[b]);
  aVoxel.adjoiningCoords[aVoxel.adjoiningSize++] = b;
  adjoin.adjoiningCoords[adjoin.adjoiningSize++] = a;
}

void MicrotubuleProcess::connectSeamEastWest(unsigned i)
{
  unsigned a(startCoord+i);
  Voxel& aVoxel((*theLattice)[a]);
  unsigned b(startCoord+(Protofilaments-1)*theDimerSize+i); 
  Voxel& adjoin((*theLattice)[b]);
  aVoxel.adjoiningCoords[aVoxel.adjoiningSize++] = b;
  adjoin.adjoiningCoords[adjoin.adjoiningSize++] = a;
}

void MicrotubuleProcess::connectNwSw(unsigned i)
{
  unsigned a(startCoord+i);
  Voxel& aVoxel((*theLattice)[a]);
  unsigned b(startCoord+(Protofilaments-1)*theDimerSize+(i-1)); 
  Voxel& adjoin((*theLattice)[b]);
  aVoxel.adjoiningCoords[aVoxel.adjoiningSize++] = b;
  adjoin.adjoiningCoords[adjoin.adjoiningSize++] = a;
}

void MicrotubuleProcess::enlistLatticeVoxels()
{
  for(unsigned n(startCoord); n != endCoord; ++n)
    {
      Voxel& offVoxel((*theLattice)[n]);
      offVoxel.diffuseSize = offVoxel.adjoiningSize;
      double rA(theSpatiocyteStepper->getMinLatticeSpace());
      if(rA < offLatticeRadius)
        {
         rA = offLatticeRadius;
        } 
      Point center(*offVoxel.point);
      unsigned aCoord(theSpatiocyteStepper->point2coord(center));
      Point cl(theSpatiocyteStepper->coord2point(aCoord));
      //theSpecies[3]->addMolecule(aVoxel.coord);
      Point bottomLeft(*offVoxel.point);
      Point topRight(*offVoxel.point);
      bottomLeft.x -= rA+center.x-cl.x+theSpatiocyteStepper->getColLength();
      bottomLeft.y -= rA+center.y-cl.y+theSpatiocyteStepper->getLayerLength();
      bottomLeft.z -= rA+center.z-cl.z+theSpatiocyteStepper->getRowLength();
      topRight.x += rA+cl.x-center.x+theSpatiocyteStepper->getColLength()*1.5;
      topRight.y += rA+cl.y-center.y+theSpatiocyteStepper->getLayerLength()*1.5;
      topRight.z += rA+cl.z-center.z+theSpatiocyteStepper->getRowLength()*1.5;
      unsigned blRow(0);
      unsigned blLayer(0);
      unsigned blCol(0);
      theSpatiocyteStepper->point2global(bottomLeft, blRow, blLayer, blCol);
      unsigned trRow(0);
      unsigned trLayer(0);
      unsigned trCol(0);
      theSpatiocyteStepper->point2global(topRight, trRow, trLayer, trCol);
      std::vector<unsigned> checkedAdjoins;
      for(unsigned i(blRow); i < trRow; ++i)
        {
          for(unsigned j(blLayer); j < trLayer; ++j)
            {
              for(unsigned k(blCol); k < trCol; ++k)
                {
                  unsigned lat(theSpatiocyteStepper->global2coord(i, j, k));
                  Voxel& latVoxel((*theLattice)[lat]);
                  if(latVoxel.id != theSpatiocyteStepper->getNullID())
                    {
                      //theSpecies[3]->addMolecule(latVoxel);
                      Point aPoint(theSpatiocyteStepper->coord2point(lat));
                      if(inMTCylinder(aPoint))
                        {
                          for(unsigned l(0); l != theAdjoiningCoordSize;
                              ++l)
                            {
                              unsigned adj(latVoxel.adjoiningCoords[l]);
                              Voxel& adjoin((*theLattice)[adj]);
                              if(adjoin.id == theComp->vacantSpecies->getID())
                                {
                                  checkedAdjoins.push_back(adj);
                                  addDirect(offVoxel, n, adjoin, adj);
                                }
                            }
                        }
                    }
                }
            }
        }
      for(unsigned i(0); i != checkedAdjoins.size(); ++i)
        {
          (*theLattice)[checkedAdjoins[i]].id = theComp->vacantSpecies->getID();
        }
    }
  for(unsigned i(0); i != occCoords.size(); ++i)
    {
      Voxel& aVoxel((*theLattice)[occCoords[i]]);
      unsigned* temp = aVoxel.initAdjoins;
      aVoxel.initAdjoins = aVoxel.adjoiningCoords;
      aVoxel.adjoiningCoords = temp;
      aVoxel.diffuseSize = aVoxel.adjoiningSize;
    }
  for(unsigned i(0); i != occCoords.size(); ++i)
    {
      Voxel& aVoxel((*theLattice)[occCoords[i]]);
      for(unsigned i(0); i != aVoxel.adjoiningSize; ++i)
        {
          unsigned aCoord(aVoxel.adjoiningCoords[i]);
          Voxel& adjoin((*theLattice)[aCoord]);
          if(adjoin.id == theComp->vacantSpecies->getID())
            {
              Point adPoint(theSpatiocyteStepper->coord2point(aCoord));
              if(inMTCylinder(adPoint))
                {
                  std::cout << "error in MT Process" << std::endl;
                }
            }
          else if(adjoin.id != theVacantSpecies->getID() &&
                  adjoin.id != theMinusSpecies->getID() &&
                  adjoin.id != thePlusSpecies->getID() &&
                  adjoin.id != theSpatiocyteStepper->getNullID())
            {
              std::cout << "species error in MT Process" << std::endl;
            }
        }
    }
}

void MicrotubuleProcess::addDirect(Voxel& offVoxel, unsigned a,
                                   Voxel& adjoin, unsigned b)
{
  Point aPoint(*offVoxel.point);
  adjoin.id = tempID;
  Point adPoint(theSpatiocyteStepper->coord2point(b));
  if(!inMTCylinder(adPoint))
    { 
      if(initAdjoins(adjoin)) 
        {
          occCoords.push_back(b);
        }
      double dist(getDistance(&aPoint, &adPoint)); 
      if(dist <= latticeRadius+offLatticeRadius)
        {
          offVoxel.adjoiningCoords[offVoxel.adjoiningSize++] = b;
          updateAdjoinSize(adjoin);
          adjoin.initAdjoins[adjoin.adjoiningSize++] = a;
        }
      else
        { 
          addIndirect(offVoxel, a, adjoin, b);
        }
    }
}

void MicrotubuleProcess::addIndirect(Voxel& offVoxel, unsigned a,
                                     Voxel& latVoxel, unsigned b)
{
  Point aPoint(*offVoxel.point);
  for(unsigned i(0); i != theAdjoiningCoordSize; ++i)
    {
      unsigned aCoord(latVoxel.adjoiningCoords[i]);
      Voxel& adjoin((*theLattice)[aCoord]);
      if(adjoin.id == theComp->vacantSpecies->getID() || adjoin.id == tempID)
        {
          Point adPoint(theSpatiocyteStepper->coord2point(aCoord));
          double dist(getDistance(&aPoint, &adPoint)); 
          if(dist <= offLatticeRadius && inMTCylinder(adPoint))
            { 
              offVoxel.adjoiningCoords[offVoxel.adjoiningSize++] = b; 
              initAdjoins(latVoxel);
              updateAdjoinSize(latVoxel);
              latVoxel.initAdjoins[latVoxel.adjoiningSize++] = a;
            }
        }
    }
}

bool MicrotubuleProcess::initAdjoins(Voxel& aVoxel)
{
  if(aVoxel.initAdjoins == NULL)
    {
      aVoxel.adjoiningSize = 0;
      aVoxel.initAdjoins = new unsigned[theAdjoiningCoordSize];
      for(unsigned i(0); i != theAdjoiningCoordSize; ++i)
        {
          unsigned aCoord(aVoxel.adjoiningCoords[i]);
          Point aPoint(theSpatiocyteStepper->coord2point(aCoord));
          if(!inMTCylinder(aPoint))
            {
              aVoxel.initAdjoins[aVoxel.adjoiningSize++] = aCoord;
            }
        }
      return true;
    }
  return false;
}

void MicrotubuleProcess::updateAdjoinSize(Voxel& aVoxel)
{
 if(aVoxel.adjoiningSize >= theAdjoiningCoordSize)
    {
      unsigned* temp(new unsigned[aVoxel.adjoiningSize+1]);
      for(unsigned i(0); i != aVoxel.adjoiningSize; ++i)
        {
          temp[i] = aVoxel.initAdjoins[i];
        }
      delete[] aVoxel.initAdjoins;
      aVoxel.initAdjoins = temp;
    }
}


bool MicrotubuleProcess::inMTCylinder(Point& N)
{
  Point E(M);
  Point W(P);
  Point S(M);
  double t((-E.x*N.x-E.y*N.y-E.z*N.z+E.x*S.x+E.y*S.y+E.z*S.z+N.x*W.x-S.x*W.x+N.y*W.y-S.y*W.y+N.z*W.z-S.z*W.z)/(E.x*E.x+E.y*E.y+E.z*E.z-2*E.x*W.x+W.x*W.x-2*E.y*W.y+W.y*W.y-2*E.z*W.z+W.z*W.z));
  double dist(sqrt(pow(-N.x+S.x+t*(-E.x+W.x),2)+pow(-N.y+S.y+t*(-E.y+W.y),2)+pow(-N.z+S.z+t*(-E.z+W.z),2)));
  if(dist < nSubunitRadius)
    {
      return true;
    }
  return false;
}
*/





