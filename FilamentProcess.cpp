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

#include "FilamentProcess.hpp"
#include "Vector.hpp"

LIBECS_DM_INIT(FilamentProcess, Process); 

unsigned FilamentProcess::getLatticeResizeCoord(unsigned aStartCoord)
{
  theComp = theSpatiocyteStepper->system2Comp(getSuperSystem());
  theVacantSpecies->resetFixedAdjoins();
  theVacantSpecies->setRadius(SubunitRadius);
  tempID = theSpecies.size();
  //The compartment center point (origin):
  C = theComp->centerPoint;
  C.x += OriginX*theComp->lengthX/2;
  C.y += OriginY*theComp->lengthY/2;
  C.z += OriginZ*theComp->lengthZ/2;
  setCompartmentDimension();
  for(unsigned i(0); i != theFilamentSpecies.size(); ++i)
    {
      theFilamentSpecies[i]->setIsOffLattice();
      theFilamentSpecies[i]->setDimension(dimension);
      theFilamentSpecies[i]->setVacantSpecies(theVacantSpecies);
      theFilamentSpecies[i]->setRadius(SubunitRadius);
    }
  //Normalized compartment lengths in terms of lattice voxel radius:
  normLength = Length/(VoxelRadius*2);
  normWidth = Width/(VoxelRadius*2);
  startCoord = aStartCoord;
  endCoord = startCoord+Filaments*Subunits;
  return endCoord-startCoord;
}

void FilamentProcess::setCompartmentDimension()
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

void FilamentProcess::initializeThird()
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

void FilamentProcess::initializeFilaments()
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
  Point D; //The direction vector from the minus end to the random point, R
  D.x = R.x-M.x;
  D.y = R.y-M.y;
  D.z = R.z-M.z;
  double NormD(sqrt(D.x*D.x+D.y*D.y+D.z*D.z));
  D.x /= NormD;
  D.y /= NormD;
  D.z /= NormD;
  //std::cout << "D.x:" << D.x << " y:" << D.y << " z:" << D.z << std::endl;
  //std::cout << "T.x:" << T.x << " y:" << T.y << " z:" << T.z << std::endl;
  //The start point of the first protofilament:
  Point S; 
  S.x = M.x;
  S.y = M.y;
  S.z = M.z;
  //std::cout << "S.x:" << S.x << " y:" << S.y << " z:" << S.z << std::endl;
  addCompVoxel(0, 0, S);
  for(unsigned i(1); i != Filaments; ++i)
    {
      Point U(S);
      U.x += i*normSubunitRadius*sqrt(3)*D.x;
      U.y += i*normSubunitRadius*sqrt(3)*D.y;
      U.z += i*normSubunitRadius*sqrt(3)*D.z;
      if(i%2 == 1)
        {
          U.x += normSubunitRadius*T.x;
          U.y += normSubunitRadius*T.y;
          U.z += normSubunitRadius*T.z;
        }
      addCompVoxel(i, 0, U);
    }
}

void FilamentProcess::initializeDirectionVector()
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
  M.x += C.x;
  M.y += C.y;
  M.z += C.z;
  //Direction vector from the Minus end to center
  T.x = C.x-M.x;
  T.y = C.y-M.y;
  T.z = C.z-M.z;
  //Make T a unit vector
  double NormT(sqrt(T.x*T.x+T.y*T.y+T.z*T.z));
  T.x /= NormT;
  T.y /= NormT;
  T.z /= NormT;
  //Rotated Plus end
  P.x = M.x+Length*T.x;
  P.y = M.y+Length*T.y;
  P.z = M.z+Length*T.z;
}

void FilamentProcess::addCompVoxel(unsigned filamentIndex, 
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

void FilamentProcess::elongateFilaments()
{
  for(unsigned i(0); i != Filaments; ++i)
    {
      Voxel* startVoxel(&(*theLattice)[startCoord+i*Subunits]);
      Point A(*startVoxel->point);
      for(unsigned j(1); j != Subunits; ++j)
        {
          A.x += (normSubunitRadius*2)*T.x;
          A.y += (normSubunitRadius*2)*T.y;
          A.z += (normSubunitRadius*2)*T.z;
          addCompVoxel(i, j, A);
        }
    }
}

void FilamentProcess::connectFilaments()
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

void FilamentProcess::connectPeriodic(unsigned j)
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

void FilamentProcess::connectNorthSouth(unsigned i, unsigned j)
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

void FilamentProcess::connectEastWest(unsigned i, unsigned j)
{
  unsigned a(startCoord+j*Subunits+i);
  Voxel& aVoxel((*theLattice)[a]);
  unsigned b(startCoord+(j-1)*Subunits+i); 
  Voxel& adjoin((*theLattice)[b]);
  aVoxel.adjoiningCoords[aVoxel.adjoiningSize++] = b;
  adjoin.adjoiningCoords[adjoin.adjoiningSize++] = a;
}

void FilamentProcess::connectSeamEastWest(unsigned i)
{
  unsigned a(startCoord+i);
  Voxel& aVoxel((*theLattice)[a]);
  unsigned b(startCoord+(Filaments-1)*Subunits+i); 
  Voxel& adjoin((*theLattice)[b]);
  aVoxel.adjoiningCoords[aVoxel.adjoiningSize++] = b;
  adjoin.adjoiningCoords[adjoin.adjoiningSize++] = a;
}

void FilamentProcess::connectNwSw(unsigned i)
{
  unsigned a(startCoord+i);
  Voxel& aVoxel((*theLattice)[a]);
  unsigned b(startCoord+(Filaments-1)*Subunits+(i-1)); 
  Voxel& adjoin((*theLattice)[b]);
  aVoxel.adjoiningCoords[aVoxel.adjoiningSize++] = b;
  adjoin.adjoiningCoords[adjoin.adjoiningSize++] = a;
}

void FilamentProcess::setCompartmentVectors()
{
  Point FilamentStart(*(*theLattice)[startCoord].point);
  //FilamentStart = sub(FilamentStart, C);
  Point FilamentEnd(*(*theLattice)[startCoord+Subunits-1].point);
  //FilamentEnd = sub(FilamentEnd, C);
  subunitVector = sub(FilamentEnd, FilamentStart);
  if(dimension == 2)
    {
      Point SurfaceEnd(*(*theLattice)[startCoord+(Filaments-1)*Subunits].point);
      //SurfaceEnd = sub(SurfaceEnd, C);
      filamentVector = sub(SurfaceEnd, FilamentStart);
      surfaceNormal = cross(subunitVector, filamentVector);
      surfaceDisplacement = dot(surfaceNormal, SurfaceEnd);
    }
}

void FilamentProcess::interfaceSubunits()
{
  enlistInterfaceVoxels();
  enlistNonIntersectInterfaceVoxels();
  theSpecies[5]->setIsPopulated();
  theSpecies[6]->setIsPopulated();
  std::cout << "interface size:" << theInterfaceSpecies->size() << std::endl;
}

void FilamentProcess::enlistInterfaceVoxels()
{
  subunitInterfaces.resize(Filaments*Subunits);
  for(unsigned i(startCoord); i != endCoord; ++i)
    {
      Voxel& subunit((*theLattice)[i]);
      //theSpecies[1]->addMolecule(&subunit);
      subunit.diffuseSize = subunit.adjoiningSize;
      Point center(*subunit.point);
      Point bottomLeft(*subunit.point);
      Point topRight(*subunit.point);
      bottomLeft.x -= normSubunitRadius+theSpatiocyteStepper->getColLength();
      bottomLeft.y -= normSubunitRadius+theSpatiocyteStepper->getLayerLength();
      bottomLeft.z -= normSubunitRadius+theSpatiocyteStepper->getRowLength();
      topRight.x += normSubunitRadius+theSpatiocyteStepper->getColLength();
      topRight.y += normSubunitRadius+theSpatiocyteStepper->getLayerLength();
      topRight.z += normSubunitRadius+theSpatiocyteStepper->getRowLength();
      unsigned blRow(0);
      unsigned blLayer(0);
      unsigned blCol(0);
      theSpatiocyteStepper->point2global(bottomLeft, blRow, blLayer, blCol);
      //unsigned aCoord(theSpatiocyteStepper->point2coord(bottomLeft));
      //theSpecies[5]->addMolecule(&(*theLattice)[aCoord]);
      //unsigned bCoord(theSpatiocyteStepper->point2coord(topRight));
      //theSpecies[5]->addMolecule(&(*theLattice)[bCoord]);
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

void FilamentProcess::addInterfaceVoxel(unsigned subunitCoord,
                                        unsigned voxelCoord)
{ 
  Voxel& subunit((*theLattice)[subunitCoord]);
  Point subunitPoint(*subunit.point);
  Point voxelPoint(theSpatiocyteStepper->coord2point(voxelCoord));
  double dist(getDistance(&subunitPoint, &voxelPoint));
  if(dist <= normSubunitRadius+normVoxelRadius) 
    {
      Voxel& voxel((*theLattice)[voxelCoord]);
      //theSpecies[6]->addMolecule(&voxel);
      //Insert voxel in the list of interface voxels if was not already:
      //if(voxel.id == theComp->vacantSpecies->getID())
      if(voxel.id != theInterfaceSpecies->getID())
        {
          theSpecies[5]->addMolecule(&voxel);
          theInterfaceSpecies->addMolecule(&voxel);
        }
      subunitInterfaces[subunitCoord-startCoord].push_back(voxelCoord);
    }
}

void FilamentProcess::enlistNonIntersectInterfaceVoxels()
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
              addNonIntersectInterfaceVoxel(adjoin);
            }
        }
    }
}

void FilamentProcess::addNonIntersectInterfaceVoxel(Voxel& aVoxel)
{
  Point pointA(theSpatiocyteStepper->coord2point(aVoxel.coord));
  double distA(point2planeDist(pointA, surfaceNormal, surfaceDisplacement));
  for(unsigned i(0); i != theAdjoiningCoordSize; ++i)
    {
      Voxel& adjoin((*theLattice)[aVoxel.adjoiningCoords[i]]);
      if(adjoin.id != theInterfaceSpecies->getID())
        {
          Point pointB(theSpatiocyteStepper->coord2point(adjoin.coord));
          double distB(point2planeDist(pointB, surfaceNormal,
                                       surfaceDisplacement));
          if((distA < 0) != (distB < 0))
            {
              if(abs(distA) < abs(distB))
                { 
                  theSpecies[5]->addMolecule(&aVoxel);
                  theInterfaceSpecies->addMolecule(&aVoxel);
                  return;
                }
              else
                {
                  theSpecies[5]->addMolecule(&adjoin);
                  theInterfaceSpecies->addMolecule(&adjoin);
                }
            }
        }
    }
}

void FilamentProcess::addDirect(Voxel& subunit, unsigned a,
                                   Voxel& adjoin, unsigned b)
{
  Point aPoint(*subunit.point);
  adjoin.id = tempID;
  Point adPoint(theSpatiocyteStepper->coord2point(b));
  if(!inMTCylinder(adPoint))
    { 
      if(initAdjoins(adjoin)) 
        {
          occCoords.push_back(b);
        }
      double dist(getDistance(&aPoint, &adPoint)); 
      if(dist <= normVoxelRadius+normSubunitRadius)
        {
          subunit.adjoiningCoords[subunit.adjoiningSize++] = b;
          updateAdjoinSize(adjoin);
          adjoin.initAdjoins[adjoin.adjoiningSize++] = a;
        }
      else
        { 
          addIndirect(subunit, a, adjoin, b);
        }
    }
}

void FilamentProcess::addIndirect(Voxel& subunit, unsigned a,
                                     Voxel& latVoxel, unsigned b)
{
  Point aPoint(*subunit.point);
  for(unsigned i(0); i != theAdjoiningCoordSize; ++i)
    {
      unsigned aCoord(latVoxel.adjoiningCoords[i]);
      Voxel& adjoin((*theLattice)[aCoord]);
      if(adjoin.id == theComp->vacantSpecies->getID() || adjoin.id == tempID)
        {
          Point adPoint(theSpatiocyteStepper->coord2point(aCoord));
          double dist(getDistance(&aPoint, &adPoint)); 
          if(dist <= normSubunitRadius && inMTCylinder(adPoint))
            { 
              subunit.adjoiningCoords[subunit.adjoiningSize++] = b; 
              initAdjoins(latVoxel);
              updateAdjoinSize(latVoxel);
              latVoxel.initAdjoins[latVoxel.adjoiningSize++] = a;
            }
        }
    }
}

bool FilamentProcess::initAdjoins(Voxel& aVoxel)
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

void FilamentProcess::updateAdjoinSize(Voxel& aVoxel)
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


bool FilamentProcess::inMTCylinder(Point& N)
{
  /*
  Point E(M);
  Point W(P);
  Point S(M);
  double t((-E.x*N.x-E.y*N.y-E.z*N.z+E.x*S.x+E.y*S.y+E.z*S.z+N.x*W.x-S.x*W.x+N.y*W.y-S.y*W.y+N.z*W.z-S.z*W.z)/(E.x*E.x+E.y*E.y+E.z*E.z-2*E.x*W.x+W.x*W.x-2*E.y*W.y+W.y*W.y-2*E.z*W.z+W.z*W.z));
  double dist(sqrt(pow(-N.x+S.x+t*(-E.x+W.x),2)+pow(-N.y+S.y+t*(-E.y+W.y),2)+pow(-N.z+S.z+t*(-E.z+W.z),2)));
  if(dist < Radius)
    {
      return true;
    }
  */
  return false;
}


// The function returns the result when the point (x,y,z) is rotated about
// the line through (a,b,c) with unit direction vector ⟨u,v,w⟩ by the angle.
void FilamentProcess::rotatePointAlongVector(Point& S, double angle)
{
  double x(S.x);
  double y(S.y);
  double z(S.z);
  double a(M.x);
  double b(M.y);
  double c(M.z);
  double u(T.x);
  double v(T.y);
  double w(T.z);
  double u2(u*u);
  double v2(v*v);
  double w2(w*w);
  double cosT(cos(angle));
  double oneMinusCosT(1-cosT);
  double sinT(sin(angle));
  double xx((a*(v2 + w2) - u*(b*v + c*w - u*x - v*y - w*z)) * oneMinusCosT
                + x*cosT + (-c*v + b*w - w*y + v*z)*sinT);
  double yy((b*(u2 + w2) - v*(a*u + c*w - u*x - v*y - w*z)) * oneMinusCosT
                + y*cosT + (c*u - a*w + w*x - u*z)*sinT);
  double zz((c*(u2 + v2) - w*(a*u + b*v - u*x - v*y - w*z)) * oneMinusCosT
                + z*cosT + (-b*u + a*v - v*x + u*y)*sinT);
  S.x = xx;
  S.y = yy;
  S.z = zz;
}


/*
void FilamentProcess::enlistLatticeVoxels()
{
  interfaceVoxels.resize(Filaments*Subunits);
  directVoxels.resize(Filaments*Subunits);
  for(unsigned n(startCoord); n != endCoord; ++n)
    {
      if(n != startCoord+(Subunits*5)+2)
        {
          continue;
        }
      std::cout << "n:" << n << " subs:" << Subunits << " Fils:" << Filaments << std::endl;
      Voxel& subunit((*theLattice)[n]);
      theSpecies[1]->addMolecule(&subunit);
      subunit.diffuseSize = subunit.adjoiningSize;
      Point center(*subunit.point);
      Point bottomLeft(*subunit.point);
      Point topRight(*subunit.point);
      bottomLeft.x -= normSubunitRadius+theSpatiocyteStepper->getColLength();
      bottomLeft.y -= normSubunitRadius+theSpatiocyteStepper->getLayerLength();
      bottomLeft.z -= normSubunitRadius+theSpatiocyteStepper->getRowLength();
      topRight.x += normSubunitRadius+theSpatiocyteStepper->getColLength();
      topRight.y += normSubunitRadius+theSpatiocyteStepper->getLayerLength();
      topRight.z += normSubunitRadius+theSpatiocyteStepper->getRowLength();
      unsigned blRow(0);
      unsigned blLayer(0);
      unsigned blCol(0);
      theSpatiocyteStepper->point2global(bottomLeft, blRow, blLayer, blCol);
      std::cout << "brow:" << blRow << " blayer:" << blLayer << " bCol:" << blCol << std::endl;
      unsigned aCoord(theSpatiocyteStepper->point2coord(bottomLeft));
      //theSpatiocyteStepper->coord2global(aCoord, blRow, blLayer, blCol);
      std::cout << "brow:" << blRow << " blayer:" << blLayer << " bCol:" << blCol << std::endl;
      if(blRow > 0)
        {
          ++blRow;
        }
      if(blCol > 0)
        {
          ++blCol;
        }
      if(blLayer > 0)
        {
          ++blLayer;
        }
      theSpecies[5]->addMolecule(&(*theLattice)[aCoord]);
      unsigned bCoord(theSpatiocyteStepper->point2coord(topRight));
      theSpecies[5]->addMolecule(&(*theLattice)[bCoord]);
      unsigned trRow(0);
      unsigned trLayer(0);
      unsigned trCol(0);
      theSpatiocyteStepper->point2global(topRight, trRow, trLayer, trCol);
      std::cout << "trow:" << trRow << " tlayer:" << trLayer << " tCol:" << trCol << std::endl;
      theSpatiocyteStepper->coord2global(bCoord, trRow, trLayer, trCol);
      std::cout << "trow:" << trRow << " tlayer:" << trLayer << " tCol:" << trCol << std::endl;
      std::vector<unsigned> checkedAdjoins;
      for(unsigned i(0); i < startCoord; ++i)
        {
          setDirectAndInterface(n, i);
        }
      for(unsigned i(blRow); i < trRow; ++i)
        {
          for(unsigned j(blLayer); j < trLayer; ++j)
            {
              for(unsigned k(blCol); k < trCol; ++k)
                {
                  unsigned lat(theSpatiocyteStepper->global2coord(i, j, k));
                  setDirectAndInterface(n, lat);
                }
            }
        }
    }
  //connectInterfaceVoxels();
  theSpecies[4]->setIsPopulated();
  theSpecies[5]->setIsPopulated();
  theSpecies[1]->setIsPopulated();
}
        */

