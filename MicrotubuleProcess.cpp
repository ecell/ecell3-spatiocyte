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
  vacantVoxels.resize(Quantity);
  theProcessSpecies[9]->setIsOffLattice();
}

void MicrotubuleProcess::initializeFourth()
{
  queueStartVoxels();
  unsigned int i(0);
  for(unsigned int j(0); i != Quantity && j != Quantity*30; ++j)
    {
      Voxel* aVoxel(startVoxels[j]);
      if(!checkStartVoxel(aVoxel))
        {
          continue;
        }
      Voxel* aParent(NULL);
      addVacantVoxel(i, aVoxel);
      //theProcessSpecies[i]->addMolecule(aVoxel);
      Point S(theSpatiocyteStepper->coord2point(aVoxel->coord));
      Voxel* aNeighbor(NULL);
      int count(28);
      do {
        --count;
        double aDist;
        aNeighbor = getNeighbor(aVoxel, S, aParent, aDist); 
        if(aNeighbor)
          {
            addVacantVoxel(i, aNeighbor);
            //theProcessSpecies[i]->addMolecule(aNeighbor);
            aParent = aVoxel;
            aVoxel = aNeighbor;
          }
      } while(aNeighbor);
      Point minPoint(theSpatiocyteStepper->coord2point(vacantVoxels[i].front()->coord));
      Point maxPoint(theSpatiocyteStepper->coord2point(vacantVoxels[i].back()->coord));
      std::cout << "the distance:" << getDistance(&minPoint, &maxPoint);
      std::cout << " max distance:" << maxLength << std::endl;
      if(getDistance(&minPoint, &maxPoint) < 0.9*maxLength)
        {
          removeVacantVoxels(i);
        }
      else
        {
          ++i;
        }
    }
  for(unsigned int i(0); i != vacantVoxels.size(); ++i)
    {
      for(std::vector<Voxel*>::iterator j(vacantVoxels[i].begin());
          j != vacantVoxels[i].end(); ++j)
        { 
          theProcessSpecies[i]->addMolecule(*j);
        }
    }
}

void MicrotubuleProcess::addVacantVoxel(unsigned int anIndex, Voxel* aVoxel)
{
  vacantVoxels[anIndex].push_back(aVoxel);
  aVoxel->id = theProcessSpecies[anIndex]->getID();
}

void MicrotubuleProcess::removeVacantVoxels(unsigned int anIndex)
{
  for(std::vector<Voxel*>::iterator i(vacantVoxels[anIndex].begin());
      i != vacantVoxels[anIndex].end(); ++i)
    { 
      (*i)->id = theComp->vacantID;
    }
  vacantVoxels[anIndex].resize(0);
}

bool MicrotubuleProcess::checkStartVoxel(Voxel* aVoxel)
{
  aVoxel->id = theProcessSpecies[0]->getID();
  Point S(theSpatiocyteStepper->coord2point(aVoxel->coord));
  double aDist;
  Voxel* firstNeighbor(getNeighbor(aVoxel, S, NULL, aDist)); 
  //std::cout << "firstDist:" << aDist << std::endl;
  if(!firstNeighbor || aDist > 0.7)
    {
      return false;
    }
  else
    {
      firstNeighbor->id = theProcessSpecies[0]->getID();
    }
  Voxel* secondNeighbor(getNeighbor(firstNeighbor, S, aVoxel, aDist)); 
  //std::cout << "secondDist:" << aDist << std::endl;
  if(!secondNeighbor || aDist > 0.7)
    {
      return false;
    }
  aVoxel->id = theComp->vacantID;
  firstNeighbor->id = theComp->vacantID;
  return true;
}

Voxel* MicrotubuleProcess::getNeighbor(Voxel* aVoxel, Point& S, Voxel* aParent,
                                           double& shortestDist)
{
  shortestDist = 1e+10;
  Voxel* aNeighbor(NULL);
  Point currS(S);
  for(int j(0); j != aVoxel->adjoiningSize; ++j)
    {
      Voxel* anAdjoin(aVoxel->adjoiningVoxels[j]);
      Point N(theSpatiocyteStepper->coord2point(anAdjoin->coord));
      if(anAdjoin->id == theComp->vacantID && notNeighbor(aParent, anAdjoin))
        {
          double t((-E.x*N.x-E.y*N.y-E.z*N.z+E.x*S.x+E.y*S.y+E.z*S.z+N.x*W.x-S.x*W.x+N.y*W.y-S.y*W.y+N.z*W.z-S.z*W.z)/(E.x*E.x+E.y*E.y+E.z*E.z-2*E.x*W.x+W.x*W.x-2*E.y*W.y+W.y*W.y-2*E.z*W.z+W.z*W.z));
          if(t < -0.001)
            {
              double dist(sqrt(pow(-N.x+S.x+t*(-E.x+W.x),2)+pow(-N.y+S.y+t*(-E.y+W.y),2)+pow(-N.z+S.z+t*(-E.z+W.z),2)));

              //std::cout << "t:" << t << " dist:" << dist << std::endl;
              if(dist < shortestDist && dist < 1.0)
                {
                  Point tempS;
                  tempS.x = S.x+t*(-E.x+W.x);
                  tempS.y = S.y+t*(-E.y+W.y);
                  tempS.z = S.z+t*(-E.z+W.z);
                  if(notShared(anAdjoin, tempS, aVoxel))
                    { 
                      //std::cout << "selected dist:" << dist << std::endl;
                      aNeighbor = anAdjoin;
                      shortestDist = dist;
                      currS = tempS;
                    }
                }
            }
        }
    } 
  S = currS;
  return aNeighbor;
}

bool MicrotubuleProcess::notShared(Voxel* aVoxel, Point S, Voxel* aParent)
{
  double shortestDist(1e+10);
  Voxel* aNeighbor(NULL);
  int count(0);
  for(int j(0); j != aVoxel->adjoiningSize; ++j)
    {
      Voxel* anAdjoin(aVoxel->adjoiningVoxels[j]);
      Point N(theSpatiocyteStepper->coord2point(anAdjoin->coord));
      if(anAdjoin->id == theComp->vacantID)
        {
          double t((-E.x*N.x-E.y*N.y-E.z*N.z+E.x*S.x+E.y*S.y+E.z*S.z+N.x*W.x-S.x*W.x+N.y*W.y-S.y*W.y+N.z*W.z-S.z*W.z)/(E.x*E.x+E.y*E.y+E.z*E.z-2*E.x*W.x+W.x*W.x-2*E.y*W.y+W.y*W.y-2*E.z*W.z+W.z*W.z));
          if(t < 0.001)
            {
              ++count;
              double dist(sqrt(pow(-N.x+S.x+t*(-E.x+W.x),2)+pow(-N.y+S.y+t*(-E.y+W.y),2)+pow(-N.z+S.z+t*(-E.z+W.z),2)));
              if(dist < shortestDist && dist < 1.0)
                {
                  aNeighbor = anAdjoin;
                  shortestDist = dist;
                }
            }
        }
    } 
  if(!aNeighbor)
    {
      return true;
    }
  else if(notNeighbor(aParent, aNeighbor) || count > 1)
    {
      return true;
    }
  return false;
}

bool MicrotubuleProcess::notNeighbor(Voxel* aSource, Voxel* aTarget)
{
  if(!aSource)
    {
      return true;
    }
  for(int i(0); i != aSource->adjoiningSize; ++i)
    {
      if(aSource->adjoiningVoxels[i] == aTarget)
        {
          return false;
        }
    }
  return true;
}

void MicrotubuleProcess::initializeDirectionVector()
{
  /* Mathematica code:
   * East = {Ex, Ey, Ez};
   * West = {Wx, Wy, Wz};
   * Direct = West - East;
   * Start = {Sx, Sy, Sz};
   * LineIntersect = Start + t*Direct
   * NeighborPoint = {Nx, Ny, Nz};
   * Solve[(LineIntersect - NeighborPoint).Direct == 0, t]
   * Distance = Norm[LineIntersect - NeighborPoint]
   */
  Point C(theComp->centerPoint);
  //East point
  E.x = theComp->lengthX/2;
  E.y = 0;
  E.z = 0;
  //West point
  W.x = -theComp->lengthX/2;
  W.y = 0;
  W.z = 0;
  theSpatiocyteStepper->rotateX(theComp->rotateX, &E, -1);
  theSpatiocyteStepper->rotateY(theComp->rotateY, &E, -1);
  theSpatiocyteStepper->rotateZ(theComp->rotateZ, &E, -1);
  theSpatiocyteStepper->rotateX(theComp->rotateX, &W, -1);
  theSpatiocyteStepper->rotateY(theComp->rotateY, &W, -1);
  theSpatiocyteStepper->rotateZ(theComp->rotateZ, &W, -1);
  E.x += C.x;
  E.y += C.y;
  E.z += C.z;
  W.x += C.x;
  W.y += C.y;
  W.z += C.z;
  //Direction vector from west to east
  D.x = E.x - W.x;
  D.y = E.y - W.y;
  D.z = E.z - W.z;
  Voxel* aVoxel(new Voxel);
  aVoxel->point = &E;
  Voxel* bVoxel(new Voxel);
  bVoxel->point = &W;
  theProcessSpecies[9]->addMolecule(aVoxel);
  theProcessSpecies[9]->addMolecule(bVoxel);
}

double MicrotubuleProcess::getWestPlaneDist(Voxel* aVoxel)
{
  Point T(theSpatiocyteStepper->coord2point(aVoxel->coord));
  double dist(((W.x-T.x)*D.x+(W.y-T.y)*D.y+(W.y-T.y)*D.y)/sqrt(D.x*D.x+D.y*D.y+D.z*D.z));
  return sqrt(dist*dist);
}

void MicrotubuleProcess::queueStartVoxels()
{
  initializeDirectionVector();






  std::vector<double> voxelDists;
  double maxDist(-1);
  Voxel* maxVoxel(NULL);
  for(std::vector<unsigned int>::iterator i(theComp->coords.begin());
      i != theComp->coords.end(); ++i)
    { 
      Voxel* aVoxel(theSpatiocyteStepper->coord2voxel(*i));
      double aDist(getWestPlaneDist(aVoxel));
      if(aDist > maxDist)
        {
          maxDist = aDist;
          maxVoxel = aVoxel;
        }
      if(startVoxels.size() != Quantity*30)
        {
          startVoxels.push_back(aVoxel);
          voxelDists.push_back(aDist);
        }
      unsigned int j(voxelDists.size()-1);
      for(; j != 0; --j)
        {
          if(voxelDists[j-1] > aDist)
            {
              voxelDists[j] = voxelDists[j-1];
              startVoxels[j] = startVoxels[j-1];
            }
          else
            {
              break;
            }
        }
      if(j != voxelDists.size()-1)
        {
          voxelDists[j] = aDist;
          startVoxels[j] = aVoxel;
        }
    }
  Point minPoint(theSpatiocyteStepper->coord2point(startVoxels[0]->coord));
  Point maxPoint(theSpatiocyteStepper->coord2point(maxVoxel->coord));
  maxLength = getDistance(&minPoint, &maxPoint);
}




