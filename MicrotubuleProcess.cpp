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
  theVacantSpecies = theProcessSpecies[0];
  theVacantSpecies->setIsVacant();
  C = theComp->centerPoint;
  C.x += OriginX*theComp->lengthX/2;
  C.y += OriginY*theComp->lengthY/2;
  C.z += OriginZ*theComp->lengthZ/2;
  VoxelDiameter = theSpatiocyteStepper->getVoxelRadius()*2;
  DimerPitch /= VoxelDiameter;
  Length /= VoxelDiameter;
  MonomerPitch /= VoxelDiameter;
  Radius /= VoxelDiameter;
  theDimerSize = (unsigned int)rint(Length/DimerPitch);
  theLattice.resize(Protofilaments*theDimerSize);
  thePoints.resize(Protofilaments*theDimerSize);
  for(unsigned int i(0); i != theProcessSpecies.size(); ++i)
    {
      theProcessSpecies[i]->setIsOffLattice();
      if(i)
        {
          theProcessSpecies[i]->setVacantSpecies(theVacantSpecies);
        }
    }
  initProtofilaments();
  elongateProtofilaments();
  connectProtofilaments();
  theVacantSpecies->setIsPopulated();
}

void MicrotubuleProcess::addVacantVoxel(unsigned int protoIndex,
                                        unsigned int dimerIndex, Point& aPoint)
{
  Voxel& aVoxel(theLattice[protoIndex*theDimerSize+dimerIndex]);
  aVoxel.point = &thePoints[protoIndex*theDimerSize+dimerIndex];
  *aVoxel.point = aPoint;
  aVoxel.adjoiningVoxels = new Voxel*[theAdjoiningVoxelSize];
  aVoxel.adjoiningSize = 2;
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
  std::cout << "D.x:" << D.x << " y:" << D.y << " z:" << D.z << std::endl;
  std::cout << "T.x:" << T.x << " y:" << T.y << " z:" << T.z << std::endl;
  Point S; //The start point of the first protofilament
  S.x = M.x+Radius*D.x;
  S.y = M.y+Radius*D.y;
  S.z = M.z+Radius*D.z;
  std::cout << "S.x:" << S.x << " y:" << S.y << " z:" << S.z << std::endl;
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
  for(unsigned int i(0); i != Protofilaments; ++i)
    {
      for(unsigned int j(0); j != theDimerSize-1; ++j)
        { 
          Voxel& firstVoxel(theLattice[i*theDimerSize+j]);
          Voxel& secondVoxel(theLattice[i*theDimerSize+j+1]);
          firstVoxel.adjoiningVoxels[NORTH] = &secondVoxel;
          secondVoxel.adjoiningVoxels[SOUTH] = &firstVoxel;
        }
    }
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


