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
  C = theComp->centerPoint;
  C.x += OriginX*theComp->lengthX/2;
  C.y += OriginY*theComp->lengthY/2;
  C.z += OriginZ*theComp->lengthZ/2;
  VoxelDiameter = theSpatiocyteStepper->getVoxelRadius()*2;
  DimerPitch /= VoxelDiameter;
  Length /= VoxelDiameter;
  MonomerPitch /= VoxelDiameter;
  Radius /= VoxelDiameter;
  vacantVoxels.resize(Protofilaments);
  for(unsigned int i(0); i != theProcessSpecies.size(); ++i)
    {
      theProcessSpecies[i]->setIsOffLattice();
    }
}

void MicrotubuleProcess::initializeFourth()
{
  initProtofilaments();
  elongateProtofilaments();
}

void MicrotubuleProcess::addVacantVoxel(unsigned int anIndex, Voxel* aVoxel)
{

  vacantVoxels[anIndex].push_back(aVoxel);
  theProcessSpecies[anIndex+1]->addMolecule(aVoxel);
  //aVoxel->id = theProcessSpecies[anIndex]->getID();

}

void MicrotubuleProcess::removeVacantVoxels(unsigned int anIndex)
{
  for(std::vector<Voxel*>::iterator i(vacantVoxels[anIndex].begin());
      i != vacantVoxels[anIndex].end(); ++i)
    { 
      (*i)->id = theVacantSpecies->getID();
    }
  vacantVoxels[anIndex].resize(0);
}

void MicrotubuleProcess::initializeDirectionVector()
{ 
  /*
   * MEnd = {Mx, My, Mz};(*minus end*) 
   * PEnd = {Px, Py, Pz};(*plus end*)
   * MTAxis = (PEnd - MEnd)/Norm[PEnd - MEnd] (*direction vector along the MT
   * long axis*)
   */
  //Plus end
  P.x = Length/2;
  P.y = 0;
  P.z = 0;
  //Minus end
  M.x = -Length/2;
  M.y = 0;
  M.z = 0;
  theSpatiocyteStepper->rotateX(RotateX, &M, -1);
  theSpatiocyteStepper->rotateY(RotateY, &M, -1);
  theSpatiocyteStepper->rotateZ(RotateZ, &M, -1);
  theSpatiocyteStepper->rotateX(RotateX, &P, -1);
  theSpatiocyteStepper->rotateY(RotateY, &P, -1);
  theSpatiocyteStepper->rotateZ(RotateZ, &P, -1);
  P.x += C.x;
  P.y += C.y;
  P.z += C.z;
  M.x += C.x;
  M.y += C.y;
  M.z += C.z;
  //Direction vector from Minus to Plus end
  T.x = P.x-M.x;
  T.y = P.y-M.y;
  T.z = P.z-M.z;
  //Make T a unit vector
  double NormT(sqrt(T.x*T.x+T.y*T.y+T.z*T.z));
  T.x /= NormT;
  T.y /= NormT;
  T.z /= NormT;
  Voxel* aVoxel(new Voxel);
  aVoxel->point = &M;
  Voxel* bVoxel(new Voxel);
  bVoxel->point = &P;
  std::cout << "M.x:" << M.x << " y:" << M.y << " z:" << M.z << std::endl;
  std::cout << "P.x:" << P.x << " y:" << P.y << " z:" << P.z << std::endl;
  theProcessSpecies[9]->addMolecule(aVoxel);
  theProcessSpecies[9]->addMolecule(bVoxel);
}

void MicrotubuleProcess::initProtofilaments()
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
  std::cout << "D.x:" << D.x << " y:" << D.y << " z:" << D.z << std::endl;
  std::cout << "T.x:" << T.x << " y:" << T.y << " z:" << T.z << std::endl;
  Point S; //The start point of the first protofilament
  S.x = M.x+Radius*D.x;
  S.y = M.y+Radius*D.y;
  S.z = M.z+Radius*D.z;
  std::cout << "S.x:" << S.x << " y:" << S.y << " z:" << S.z << std::endl;
  Voxel* aVoxel(new Voxel);
  aVoxel->point = new Point;
  *aVoxel->point = S;
  addVacantVoxel(0, aVoxel);
  for(int i(1); i != Protofilaments; ++i)
    {
      double angle(2*M_PI/Protofilaments);
      rotatePointAlongVector(S, angle);
      S.x = S.x+MonomerPitch/(Protofilaments-1)*T.x;
      S.y = S.y+MonomerPitch/(Protofilaments-1)*T.y;
      S.z = S.z+MonomerPitch/(Protofilaments-1)*T.z;
      Voxel* bVoxel(new Voxel);
      bVoxel->point = new Point;
      *bVoxel->point = S;
      addVacantVoxel(i, bVoxel);
    }
}

void MicrotubuleProcess::elongateProtofilaments()
{
  for(unsigned int i(0); i != Protofilaments; ++i)
    {
      Voxel* startVoxel(vacantVoxels[i][0]);
      Point A(*startVoxel->point);
      for(unsigned int j(0); j != (unsigned int)rint(Length/DimerPitch); ++j)
        {
          A.x = A.x+DimerPitch*T.x;
          A.y = A.y+DimerPitch*T.y;
          A.z = A.z+DimerPitch*T.z;
          Voxel* aVoxel(new Voxel);
          aVoxel->point = new Point;
          *aVoxel->point = A;
          addVacantVoxel(i, aVoxel);
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
  std::cout << "S.x:" << S.x << " y:" << S.y << " z:" << S.z << std::endl;
}
