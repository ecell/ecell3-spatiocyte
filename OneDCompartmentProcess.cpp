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

#include "OneDCompartmentProcess.hpp"

LIBECS_DM_INIT(OneDCompartmentProcess, Process); 

void OneDCompartmentProcess::initializeThird()
{
  theComp = theSpatiocyteStepper->system2Comp(getSuperSystem());
  //VacantVoxels[0].push_back(getBeginVoxel());
}

void OneDCompartmentProcess::initializeFourth()
{
  queueStartVoxels();
  unsigned int i(0);
  for(unsigned int j(0); i != Quantity && j != Quantity*3; ++j)
    {
      Voxel* aVoxel(startVoxels[j]);
      std::cout << std::endl << "checking j:" << j << std::endl;
      if(!checkStartVoxel(aVoxel))
        {
          continue;
        }
      std::cout << "done checking j:" << j << std::endl;
      Voxel* aParent(NULL);
      theProcessSpecies[i]->addMolecule(aVoxel);
      Point S(theSpatiocyteStepper->coord2point(aVoxel->coord));
      Voxel* aNeighbor(NULL);
      int count(5);
      do {
        std::cout << "count:" << count << std::endl;
        --count;
        double aDist;
        aNeighbor = getNeighbor(aVoxel, S, aParent, aDist); 
        if(aNeighbor)
          {
            theProcessSpecies[i]->addMolecule(aNeighbor);
            aParent = aVoxel;
            aVoxel = aNeighbor;
          }
      } while(aNeighbor);
      ++i;
    }
  theProcessSpecies[0]->setIsPopulated();
}

bool OneDCompartmentProcess::checkStartVoxel(Voxel* aVoxel)
{
  aVoxel->id = theProcessSpecies[0]->getID();
  Point S(theSpatiocyteStepper->coord2point(aVoxel->coord));
  double aDist;
  Voxel* firstNeighbor(getNeighbor(aVoxel, S, NULL, aDist)); 
  std::cout << "firstDist:" << aDist << std::endl;
  if(!firstNeighbor || aDist > 0.7)
    {
      return false;
    }
  else
    {
      firstNeighbor->id = theProcessSpecies[0]->getID();
    }
  Voxel* secondNeighbor(getNeighbor(firstNeighbor, S, aVoxel, aDist)); 
  std::cout << "secondDist:" << aDist << std::endl;
  if(!secondNeighbor || aDist > 0.7)
    {
      return false;
    }
  aVoxel->id = theComp->vacantID;
  firstNeighbor->id = theComp->vacantID;
  return true;
}

Voxel* OneDCompartmentProcess::getNeighbor(Voxel* aVoxel, Point& S, Voxel* aParent,
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
          if(t < 0)
            {
              double dist(sqrt(pow(-N.x+S.x+t*(-E.x+W.x),2)+pow(-N.y+S.y+t*(-E.y+W.y),2)+pow(-N.z+S.z+t*(-E.z+W.z),2)));
              if(dist < shortestDist)
                {
                  Point tempS;
                  tempS.x = S.x+t*(-E.x+W.x);
                  tempS.y = S.y+t*(-E.y+W.y);
                  tempS.z = S.z+t*(-E.z+W.z);
                  if(notShared(anAdjoin, tempS, aVoxel))
                    { 
                      std::cout << "dist:" << dist << std::endl;
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

bool OneDCompartmentProcess::notShared(Voxel* aVoxel, Point S, Voxel* aParent)
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
          if(t < 0)
            {
              ++count;
              double dist(sqrt(pow(-N.x+S.x+t*(-E.x+W.x),2)+pow(-N.y+S.y+t*(-E.y+W.y),2)+pow(-N.z+S.z+t*(-E.z+W.z),2)));
              if(dist < shortestDist)
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

bool OneDCompartmentProcess::notNeighbor(Voxel* aSource, Voxel* aTarget)
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

void OneDCompartmentProcess::initializeDirectionVector()
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
}

double OneDCompartmentProcess::getWestPlaneDist(Voxel* aVoxel)
{
  Point T(theSpatiocyteStepper->coord2point(aVoxel->coord));
  double dist(((W.x-T.x)*D.x+(W.y-T.y)*D.y+(W.y-T.y)*D.y)/sqrt(D.x*D.x+D.y*D.y+D.z*D.z));
  return sqrt(dist*dist);
}

void OneDCompartmentProcess::queueStartVoxels()
{
  if(Quantity*3 > theComp->coords.size())
    {
      THROW_EXCEPTION(ValueError, String(getPropertyInterface().getClassName()) +
                                  "[" + getFullID().asString() + 
                                  "]: Quantity is larger than available compartment " +
                                  "voxels.");
    }
  initializeDirectionVector();
  std::vector<double> voxelDists;
  for(std::vector<unsigned int>::iterator i(theComp->coords.begin());
      i != theComp->coords.end(); ++i)
    { 
      Voxel* aVoxel(theSpatiocyteStepper->coord2voxel(*i));
      double aDist(getWestPlaneDist(aVoxel));
      if(startVoxels.size() != Quantity*3)
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
}




