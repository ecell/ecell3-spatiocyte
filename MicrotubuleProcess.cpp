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

void MicrotubuleProcess::initializeThird()
{
  theComp = theSpatiocyteStepper->system2Comp(getSuperSystem());
  theVacantSpecies->setIsOffLattice();
  theVacantSpecies->setRadius(DimerPitch/2);
  tempID = theSpecies.size();
  C = theComp->centerPoint;
  C.x += OriginX*theComp->lengthX/2;
  C.y += OriginY*theComp->lengthY/2;
  C.z += OriginZ*theComp->lengthZ/2;
  for(unsigned int i(0); i != theKinesinSpecies.size(); ++i)
    {
      theKinesinSpecies[i]->setIsOffLattice();
      theKinesinSpecies[i]->setVacantSpecies(theVacantSpecies);
      //theKinesinSpecies[i]->setRadius(DimerPitch/2);
    }
  VoxelDiameter = theSpatiocyteStepper->getVoxelRadius()*2;
  DimerPitch /= VoxelDiameter;
  Length /= VoxelDiameter;
  MonomerPitch /= VoxelDiameter;
  Radius /= VoxelDiameter;
  offLatticeRadius = DimerPitch/2;
  latticeRadius = 0.5;
  theDimerSize = (unsigned int)rint(Length/DimerPitch);
  theLattice.resize(Protofilaments*theDimerSize);
  thePoints.resize(Protofilaments*theDimerSize);
  initProtofilaments();
  elongateProtofilaments();
  connectProtofilaments();
  theVacantSpecies->setIsPopulated();
  connectLatticeVoxels();
}

void MicrotubuleProcess::addVacantVoxel(unsigned int protoIndex,
                                        unsigned int dimerIndex, Point& aPoint)
{
  Voxel& aVoxel(theLattice[protoIndex*theDimerSize+dimerIndex]);
  aVoxel.point = &thePoints[protoIndex*theDimerSize+dimerIndex];
  *aVoxel.point = aPoint;
  aVoxel.adjoiningVoxels = new Voxel*[theAdjoiningVoxelSize];
  aVoxel.diffuseSize = 2;
  /*
  aVoxel.adjoiningSize = 2;
  */
  for(unsigned int i(0); i != theAdjoiningVoxelSize; ++i)
    {
      aVoxel.adjoiningVoxels[i] = theNullVoxel;
    }
  theVacantSpecies->hardAddMolecule(&aVoxel);
}

void MicrotubuleProcess::initializeDirectionVector()
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
  /*
  Voxel* aVoxel(new Voxel);
  aVoxel->point = &M;
  Voxel* bVoxel(new Voxel);
  bVoxel->point = &P;
  std::cout << "M.x:" << M.x << " y:" << M.y << " z:" << M.z << std::endl;
  std::cout << "P.x:" << P.x << " y:" << P.y << " z:" << P.z << std::endl;
  theVacantSpecies->hardAddMolecule(aVoxel);
  theVacantSpecies->hardAddMolecule(bVoxel);
  */
}

void MicrotubuleProcess::initProtofilaments()
{
  theAdjoiningVoxelSize = theSpatiocyteStepper->getAdjoiningVoxelSize();
  theNullVoxel = new Voxel;
  theNullVoxel->id = theSpatiocyteStepper->getNullID();
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
  Point S; //The start point of the first protofilament
  S.x = M.x+Radius*D.x;
  S.y = M.y+Radius*D.y;
  S.z = M.z+Radius*D.z;
  //std::cout << "S.x:" << S.x << " y:" << S.y << " z:" << S.z << std::endl;
  addVacantVoxel(0, 0, S);
  for(int i(1); i != Protofilaments; ++i)
    {
      double angle(2*M_PI/Protofilaments);
      rotatePointAlongVector(S, angle);
      S.x += MonomerPitch/(Protofilaments-1)*T.x;
      S.y += MonomerPitch/(Protofilaments-1)*T.y;
      S.z += MonomerPitch/(Protofilaments-1)*T.z;
      addVacantVoxel(i, 0, S);
    }
}

void MicrotubuleProcess::elongateProtofilaments()
{
  for(unsigned int i(0); i != Protofilaments; ++i)
    {
      Voxel& startVoxel(theLattice[i*theDimerSize]);
      Point A(*startVoxel.point);
      for(unsigned int j(1); j != theDimerSize; ++j)
        {
          A.x += DimerPitch*T.x;
          A.y += DimerPitch*T.y;
          A.z += DimerPitch*T.z;
          addVacantVoxel(i, j, A);
        }
    }
}

void MicrotubuleProcess::connectProtofilaments()
{
  for(unsigned int i(0); i != theDimerSize; ++i)
    {
      for(unsigned int j(0); j != Protofilaments; ++j)
        { 
          if(i > 0)
            { 
              connectNorthSouth(i, j);
            }
          /*
          if(j > 0)
            {
              connectEastWest(i, j);
            }
            */
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
        */
    }
}

void MicrotubuleProcess::connectNorthSouth(unsigned int i, unsigned int j)
{
  Voxel* aVoxel(&theLattice[j*theDimerSize+i]);
  Voxel* adjoin(&theLattice[j*theDimerSize+(i-1)]); 
  aVoxel->adjoiningVoxels[NORTH] = adjoin;
  adjoin->adjoiningVoxels[SOUTH] = aVoxel;
  aVoxel->adjoiningSize = 2;
  adjoin->adjoiningSize = 2;
}

void MicrotubuleProcess::connectEastWest(unsigned int i, unsigned int j)
{
  Voxel* aVoxel(&theLattice[j*theDimerSize+i]);
  Voxel* adjoin(&theLattice[(j-1)*theDimerSize+i]); 
  aVoxel->adjoiningVoxels[aVoxel->adjoiningSize++] = adjoin;
  adjoin->adjoiningVoxels[adjoin->adjoiningSize++] = aVoxel;
}

void MicrotubuleProcess::connectSeamEastWest(unsigned int i)
{
  Voxel* aVoxel(&theLattice[i]);
  Voxel* adjoin(&theLattice[(Protofilaments-1)*theDimerSize+i]); 
  aVoxel->adjoiningVoxels[aVoxel->adjoiningSize++] = adjoin;
  adjoin->adjoiningVoxels[adjoin->adjoiningSize++] = aVoxel;
}

void MicrotubuleProcess::connectNwSw(unsigned int i)
{
  Voxel* aVoxel(&theLattice[i]);
  Voxel* adjoin(&theLattice[(Protofilaments-1)*theDimerSize+(i-1)]); 
  aVoxel->adjoiningVoxels[aVoxel->adjoiningSize++] = adjoin;
  adjoin->adjoiningVoxels[adjoin->adjoiningSize++] = aVoxel;
}

void MicrotubuleProcess::connectLatticeVoxels()
{
  enlistLatticeVoxels();
  /*
  theSpecies[2]->setIsPopulated();
  */
  theSpecies[3]->setIsPopulated();
}

void MicrotubuleProcess::enlistLatticeVoxels()
{
  for(std::vector<Voxel>::iterator n(theLattice.begin()); 
      n != theLattice.end(); ++n)
    {
      Voxel* offVoxel(&(*n));
      offVoxel->diffuseSize = offVoxel->adjoiningSize;
      double rA(theSpatiocyteStepper->getMinLatticeSpace());
      if(rA < offLatticeRadius)
        {
         rA = offLatticeRadius;
        } 
      Point center(*offVoxel->point);
      Voxel* aVoxel(theSpatiocyteStepper->point2voxel(center));
      Point cl(theSpatiocyteStepper->coord2point(aVoxel->coord));
      //theSpecies[3]->addMolecule(aVoxel);
      Point bottomLeft(*offVoxel->point);
      Point topRight(*offVoxel->point);
      bottomLeft.x -= rA+center.x-cl.x+theSpatiocyteStepper->getColLength();
      bottomLeft.y -= rA+center.y-cl.y+theSpatiocyteStepper->getLayerLength();
      bottomLeft.z -= rA+center.z-cl.z+theSpatiocyteStepper->getRowLength();
      topRight.x += rA+cl.x-center.x+theSpatiocyteStepper->getColLength()*1.5;
      topRight.y += rA+cl.y-center.y+theSpatiocyteStepper->getLayerLength()*1.5;
      topRight.z += rA+cl.z-center.z+theSpatiocyteStepper->getRowLength()*1.5;
      unsigned int blRow(0);
      unsigned int blLayer(0);
      unsigned int blCol(0);
      theSpatiocyteStepper->point2global(bottomLeft, blRow, blLayer, blCol);
      unsigned int trRow(0);
      unsigned int trLayer(0);
      unsigned int trCol(0);
      theSpatiocyteStepper->point2global(topRight, trRow, trLayer, trCol);
      std::vector<Voxel*> checkedAdjoins;
      for(unsigned int i(blRow); i < trRow; ++i)
        {
          for(unsigned int j(blLayer); j < trLayer; ++j)
            {
              for(unsigned int k(blCol); k < trCol; ++k)
                {
                  Voxel* latVoxel(theSpatiocyteStepper->global2voxel(i, j, k));
                  if(latVoxel->id != theSpatiocyteStepper->getNullID())
                    {
                      //theSpecies[3]->addMolecule(latVoxel);
                      Point aPoint(
                           theSpatiocyteStepper->coord2point(latVoxel->coord));
                      if(inMTCylinder(aPoint))
                        {
                          for(unsigned int l(0); l != theAdjoiningVoxelSize;
                              ++l)
                            {
                              Voxel* adjoin(latVoxel->adjoiningVoxels[l]);
                              if(adjoin->id == theComp->vacantSpecies->getID())
                                {
                                  checkedAdjoins.push_back(adjoin);
                                  addDirect(offVoxel, adjoin);
                                }
                            }
                        }
                    }
                }
            }
        }
      for(unsigned int i(0); i != checkedAdjoins.size(); ++i)
        {
          checkedAdjoins[i]->id = theComp->vacantSpecies->getID();
        }
    }
  for(unsigned int i(0); i != occVoxels.size(); ++i)
    {
      Voxel* aVoxel(occVoxels[i]);
      Voxel** temp = aVoxel->initAdjoins;
      aVoxel->initAdjoins = aVoxel->adjoiningVoxels;
      aVoxel->adjoiningVoxels = temp;
      aVoxel->diffuseSize = aVoxel->adjoiningSize;
    }
  for(unsigned int i(0); i != occVoxels.size(); ++i)
    {
      Voxel* aVoxel(occVoxels[i]);
      for(unsigned int i(0); i != aVoxel->adjoiningSize; ++i)
        {
          Voxel* adjoin(aVoxel->adjoiningVoxels[i]);
          if(adjoin->id == theComp->vacantSpecies->getID())
            {
              Point adPoint(theSpatiocyteStepper->coord2point(adjoin->coord));
              if(inMTCylinder(adPoint))
                {
                  std::cout << "error in MT Process" << std::endl;
                }
            }
          else if(adjoin->id != theVacantSpecies->getID() &&
                  adjoin->id != theSpatiocyteStepper->getNullID())
            {
              std::cout << "species error in MT Process" << std::endl;
            }
        }
    }
}

void MicrotubuleProcess::addDirect(Voxel* offVoxel, Voxel* adjoin)
{
  Point aPoint(*offVoxel->point);
  adjoin->id = tempID;
  Point adPoint(theSpatiocyteStepper->coord2point(adjoin->coord));
  if(!inMTCylinder(adPoint))
    { 
      if(initAdjoins(adjoin)) 
        {
          occVoxels.push_back(adjoin);
        }
      double dist(getDistance(&aPoint, &adPoint)); 
      if(dist <= latticeRadius+offLatticeRadius)
        {
          //std::cout << "size:" << offVoxel->adjoiningSize << std::endl;
          offVoxel->adjoiningVoxels[offVoxel->adjoiningSize++] = adjoin;
          updateAdjoinSize(adjoin);
          adjoin->initAdjoins[adjoin->adjoiningSize++] = offVoxel;
        }
      else
        { 
          addIndirect(offVoxel, adjoin);
        }
    }
}

void MicrotubuleProcess::addIndirect(Voxel* offVoxel, Voxel* latVoxel)
{
  Point aPoint(*offVoxel->point);
  for(unsigned int i(0); i != theAdjoiningVoxelSize; ++i)
    {
      Voxel* adjoin(latVoxel->adjoiningVoxels[i]);
      if(adjoin->id == theComp->vacantSpecies->getID() || 
         adjoin->id == tempID)
        {
          Point adPoint(theSpatiocyteStepper->coord2point(adjoin->coord));
          double dist(getDistance(&aPoint, &adPoint)); 
          if(dist <= offLatticeRadius && inMTCylinder(adPoint))
            { 
              //std::cout << "insize:" << offVoxel->adjoiningSize << std::endl;
              offVoxel->adjoiningVoxels[offVoxel->adjoiningSize++] = latVoxel;
              initAdjoins(latVoxel);
              updateAdjoinSize(latVoxel);
              latVoxel->initAdjoins[latVoxel->adjoiningSize++] = offVoxel;
            }
        }
    }
}

bool MicrotubuleProcess::initAdjoins(Voxel* aVoxel)
{
  if(aVoxel->initAdjoins == NULL)
    {
      aVoxel->adjoiningSize = 0;
      aVoxel->initAdjoins = new Voxel*[theAdjoiningVoxelSize];
      for(unsigned int i(0); i != theAdjoiningVoxelSize; ++i)
        {
          Voxel* adjoin(aVoxel->adjoiningVoxels[i]);
          Point aPoint(theSpatiocyteStepper->coord2point(adjoin->coord));
          if(!inMTCylinder(aPoint))
            {
              aVoxel->initAdjoins[aVoxel->adjoiningSize++] = adjoin;
            }
        }
      return true;
    }
  return false;
}

void MicrotubuleProcess::updateAdjoinSize(Voxel* aVoxel)
{
 if(aVoxel->adjoiningSize >= theAdjoiningVoxelSize)
    {
      Voxel** temp = new Voxel*[aVoxel->adjoiningSize+1];
      for(unsigned int i(0); i != aVoxel->adjoiningSize; ++i)
        {
          temp[i] = aVoxel->initAdjoins[i];
        }
      delete[] aVoxel->initAdjoins;
      aVoxel->initAdjoins = temp;
    }
}


/*
void MicrotubuleProcess::connectDirect(Voxel* offVoxel, Voxel* latVoxel,
                                       Voxel* adjoin)
{
  for(unsigned int i(0); i != theAdjoiningVoxelSize; ++i)
    {
      if(adjoin->adjoiningVoxels[i] == latVoxel)
        {
          adjoin->adjoiningVoxels[i] = offVoxel;
          offVoxel->adjoiningVoxels[offVoxel->adjoiningSize++] = adjoin;
          latticeVoxels.push_back(adjoin);
          theSpecies[2]->addMolecule(adjoin);
          return;
        }
    }
  std::cout << "error in MT connectDirect" << std::endl;
}

void MicrotubuleProcess::connectInDirect(Voxel* offVoxel, Voxel* latVoxel)
{
  Point aPoint(*offVoxel->point);
  for(unsigned int i(0); i != theAdjoiningVoxelSize; ++i)
    {
      Voxel* adjoin(latVoxel->adjoiningVoxels[i]);
      if(adjoin->id == theComp->vacantSpecies->getID())
        {
          Point adPoint(theSpatiocyteStepper->coord2point(adjoin->coord));
          double dist(getDistance(&aPoint, &adPoint)); 
          if(dist <= offLatticeRadius+latticeRadius)
            { 
              latVoxel->adjoiningVoxels[i] = offVoxel;
              offVoxel->adjoiningVoxels[offVoxel->adjoiningSize++] = latVoxel;
              latticeVoxels.push_back(latVoxel);
              theSpecies[2]->addMolecule(latVoxel);
            }

        }
    }
}
*/


void MicrotubuleProcess::addAdjoinVoxels(Voxel* offVoxel, Voxel* aVoxel)
{
  Point aPoint(*offVoxel->point);
  for(unsigned int i(0); i != theAdjoiningVoxelSize; ++i)
    {
      Voxel* adjoin(aVoxel->adjoiningVoxels[i]);
      Point adPoint(theSpatiocyteStepper->coord2point(adjoin->coord));
      double aDist(getDistance(&adPoint, &aPoint));
      if(aDist > offLatticeRadius+latticeRadius && !inMTCylinder(adPoint) &&
         notInNeighbors(offVoxel, adPoint))
       {
         theSpecies[2]->addMolecule(adjoin);
       }
    }
}

bool MicrotubuleProcess::notInNeighbors(Voxel* offVoxel, Point& aPoint)
{
  for(unsigned int i(0); i != offVoxel->adjoiningSize; ++i)
    {
      Voxel* adjoin(offVoxel->adjoiningVoxels[i]);
      Point adPoint(*adjoin->point);
      double dist(getDistance(&aPoint, &adPoint));
      if(dist <= offLatticeRadius+latticeRadius)
        {
          return false;
        }
    }
  return true;
}


bool MicrotubuleProcess::addLatticeVoxel(Voxel* aVoxel, Voxel* offVoxel)
{
  /*
    if(aVoxel->id != theComp->vacantSpecies->getID() && 
     aVoxel->id != theSpecies[2]->getID())
    {
      return false;
    }
  std::cout << "true:" << std::endl;
  if(aVoxel->id == theComp->vacantSpecies->getID())
    {
      theSpecies[2]->addMolecule(aVoxel);
      latticeVoxels.push_back(aVoxel);
    }
  return true;
  */
  if(!isValidVoxel(aVoxel))
    {
      return false;
    }
  Point targetPoint(theSpatiocyteStepper->coord2point(aVoxel->coord));
  double dist(getDistance(&targetPoint, offVoxel->point));
  if(dist < offLatticeRadius)
    {
      std::cout << "less" << std::endl;
      for(unsigned int i(0); i != aVoxel->adjoiningSize; ++i) 
        {
          Voxel* anAdjoin(aVoxel->adjoiningVoxels[i]);
          if(isValidVoxel(anAdjoin) && !isInsideVoxel(anAdjoin))
            {
              for(unsigned int j(0); j != anAdjoin->adjoiningSize; ++j) 
                {
                  if(anAdjoin->adjoiningVoxels[j] == aVoxel)
                    {
                      anAdjoin->adjoiningVoxels[j] = offVoxel;
                      if(offVoxel->adjoiningSize < theAdjoiningVoxelSize)
                        {
                          offVoxel->adjoiningVoxels[offVoxel->adjoiningSize] =
                            anAdjoin;
                          ++offVoxel->adjoiningSize;
                        } 
                      theSpecies[2]->addMolecule(anAdjoin);
                    }
                }
            }
        }
      return false;
    }
  return true;
}

void MicrotubuleProcess::connectLatticeVoxel(Voxel* aVoxel, Voxel* offVoxel)
{
  removeInsideAdjoins(aVoxel);
  if(aVoxel->diffuseSize < theAdjoiningVoxelSize)
    {
      aVoxel->adjoiningVoxels[aVoxel->diffuseSize] = offVoxel;
      ++aVoxel->diffuseSize;
    }
  if(offVoxel->adjoiningSize < theAdjoiningVoxelSize)
    {
      offVoxel->adjoiningVoxels[offVoxel->adjoiningSize] = aVoxel;
      ++offVoxel->adjoiningSize;
    }
  theSpecies[2]->addMolecule(aVoxel);
}

void MicrotubuleProcess::removeInsideAdjoins(Voxel* aVoxel)
{
  std::vector<Voxel*> outsideVoxels;
  for(unsigned int i(0); i != aVoxel->adjoiningSize; ++i)
    {
      if(!isInsideVoxel(aVoxel->adjoiningVoxels[i]))
        {
          outsideVoxels.push_back(aVoxel->adjoiningVoxels[i]);
        }
    }
  for(unsigned int i(0); i != outsideVoxels.size(); ++i)
    {
      aVoxel->adjoiningVoxels[i] = outsideVoxels[i];
    }
  aVoxel->adjoiningSize = outsideVoxels.size();
  aVoxel->diffuseSize = outsideVoxels.size();
}

bool MicrotubuleProcess::isValidVoxel(Voxel* aVoxel)
{
  if(aVoxel->id != theComp->vacantSpecies->getID() && 
     aVoxel->id != theSpecies[2]->getID())
    {
      return false;
    }
  return true;
}

bool MicrotubuleProcess::inMTCylinder(Point& N)
{
  Point E(M);
  Point W(P);
  Point S(M);
  double t((-E.x*N.x-E.y*N.y-E.z*N.z+E.x*S.x+E.y*S.y+E.z*S.z+N.x*W.x-S.x*W.x+N.y*W.y-S.y*W.y+N.z*W.z-S.z*W.z)/(E.x*E.x+E.y*E.y+E.z*E.z-2*E.x*W.x+W.x*W.x-2*E.y*W.y+W.y*W.y-2*E.z*W.z+W.z*W.z));
  double dist(sqrt(pow(-N.x+S.x+t*(-E.x+W.x),2)+pow(-N.y+S.y+t*(-E.y+W.y),2)+pow(-N.z+S.z+t*(-E.z+W.z),2)));
  if(dist < Radius)
    {
      return true;
    }
  return false;
}

bool MicrotubuleProcess::isInsideVoxel(Voxel* aVoxel)
{
  Point N(theSpatiocyteStepper->coord2point(aVoxel->coord));
  Point E(M);
  Point W(P);
  Point S(M);
  double t((-E.x*N.x-E.y*N.y-E.z*N.z+E.x*S.x+E.y*S.y+E.z*S.z+N.x*W.x-S.x*W.x+N.y*W.y-S.y*W.y+N.z*W.z-S.z*W.z)/(E.x*E.x+E.y*E.y+E.z*E.z-2*E.x*W.x+W.x*W.x-2*E.y*W.y+W.y*W.y-2*E.z*W.z+W.z*W.z));
  double dist(sqrt(pow(-N.x+S.x+t*(-E.x+W.x),2)+pow(-N.y+S.y+t*(-E.y+W.y),2)+pow(-N.z+S.z+t*(-E.z+W.z),2)));
  if(dist < Radius)
    {
      return true;
    }
  return false;
}


/*
 * The function returns the result when the point (x,y,z) is rotated about the line through (a,b,c) with unit direction vector ⟨u,v,w⟩ by the angle θ.
 * */
void MicrotubuleProcess::rotatePointAlongVector(Point& S, double angle)
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




