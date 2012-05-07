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
  VacantVoxels.resize(Quantity);
  theComp = theSpatiocyteStepper->system2Comp(getSuperSystem());
  //VacantVoxels[0].push_back(getBeginVoxel());
}

void OneDCompartmentProcess::initializeFourth()
{
  for(int j(0); j != 3; ++j)
    {
      Voxel* currVoxel(getBeginVoxel());
      if(currVoxel)
        {
          theProcessSpecies[j]->addMolecule(currVoxel);
          Point S(theSpatiocyteStepper->coord2point(currVoxel->coord)); 
          for(std::vector<unsigned int>::iterator i(theComp->coords.begin());
              i != theComp->coords.end(); ++i)
            { 
              Voxel* aVoxel(theSpatiocyteStepper->coord2voxel(*i));
              if(aVoxel->id == theComp->vacantID && isInsidePlane(aVoxel, S))
                { 
                  if(isLine(aVoxel, S))
                    {
                      theProcessSpecies[j]->addMolecule(aVoxel);
                    }
                }
            }
        }
    }
}


/*

  for(int i(0); i != 1; ++i)
    {
      Voxel* currVoxel(getBeginVoxel());
      Voxel* prevVoxel;
      do
        { 
          theProcessSpecies[i]->addMolecule(currVoxel);
          Point S(theSpatiocyteStepper->coord2point(currVoxel->coord)); 
          prevVoxel = currVoxel;
          for(int j(0); j != currVoxel->adjoiningSize; ++j)
            {
              Voxel* anAdjoin(currVoxel->adjoiningVoxels[j]); 
              if(anAdjoin->id == theComp->vacantID && isInsidePlane(anAdjoin, S))
                {
                  if(isLine(anAdjoin, S))
                    {
                      theProcessSpecies[i]->addMolecule(anAdjoin);
                      currVoxel = anAdjoin;
                      break;
                    }
                }
            }
        } while(currVoxel != prevVoxel);
    }
}
*/

bool OneDCompartmentProcess::isInsidePlane(Voxel* aVoxel, Point& S)
{ 
  Point N(theSpatiocyteStepper->coord2point(aVoxel->coord));
  double distance((-N.z+S.z)*(E.y*S.x-E.x*S.y-E.y*W.x+S.y*W.x+E.x*W.y-S.x*W.y)+
                  (-N.y+S.y)*(-E.z*S.x+E.x*S.z+E.z*W.x-S.z*W.x-E.x*W.z+S.x*W.z)+
                  (-N.z+S.x)*(E.z*S.y-E.y*S.z-E.z*W.y+S.z*W.y+E.y*W.z-S.y*W.z));
  std::cout << "distPlane:" << distance << std::endl;
  if(distance <= 0)
    {
      return true;
    }
  return false;
}

bool OneDCompartmentProcess::isLine(Voxel* aVoxel, Point& S)
{
  for(int i(0); i != aVoxel->adjoiningSize; ++i)
    { 
      Voxel* anAdjoin(aVoxel->adjoiningVoxels[i]);
      Point N(theSpatiocyteStepper->coord2point(anAdjoin->coord));
      double distance((-N.z+S.z)*(E.y*S.x-E.x*S.y-E.y*W.x+S.y*W.x+E.x*W.y-S.x*W.y)+
                      (-N.y+S.y)*(-E.z*S.x+E.x*S.z+E.z*W.x-S.z*W.x-E.x*W.z+S.x*W.z)+
                      (-N.z+S.x)*(E.z*S.y-E.y*S.z-E.z*W.y+S.z*W.y+E.y*W.z-S.y*W.z));
      std::cout << "distLine:" << distance << std::endl;
      if(distance > 0)
        {
          return true;
        }
    }
  return false;
}

/*
void OneDCompartmentProcess::initializeFourth()
{
  for(int i(0); i != 6; ++i)
    {
      Voxel* aParent(NULL);
      Voxel* aVoxel(getBeginVoxel());
      theProcessSpecies[i]->addMolecule(aVoxel);
      Point S(theSpatiocyteStepper->coord2point(aVoxel->coord));
      Voxel* aNeighbor(NULL);
      int count(8);
      do {
        --count;
        aNeighbor = getNeighbor(aVoxel, S, aParent); 
        if(aNeighbor)
          {
            theProcessSpecies[i]->addMolecule(aNeighbor);
            aParent = aVoxel;
            aVoxel = aNeighbor;
          }
      } while(aNeighbor);
    }
  theProcessSpecies[0]->setIsPopulated();
}
*/

Voxel* OneDCompartmentProcess::getNeighbor(Voxel* aVoxel, Point& S, Voxel* aParent)
{
  double shortestDist(1e+10);
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

Voxel* OneDCompartmentProcess::getBeginVoxel()
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
  double shortestDist(1e+10);
  Voxel* aBeginVoxel(NULL);
  for(std::vector<unsigned int>::iterator i(theComp->coords.begin());
      i != theComp->coords.end(); ++i)
    { 
      Voxel* aVoxel(theSpatiocyteStepper->coord2voxel(*i));
      if(aVoxel->id == theComp->vacantID)
        {
          Point T(theSpatiocyteStepper->coord2point(*i +
                                      theSpatiocyteStepper->getStartCoord()));
          double dist(((W.x-T.x)*D.x+(W.y-T.y)*D.y+(W.y-T.y)*D.y)/
                      sqrt(D.x*D.x+D.y*D.y+D.z*D.z));
          dist = sqrt(dist*dist);
          if(dist < shortestDist)
            {
              aBeginVoxel = aVoxel;
              shortestDist = dist;
            }
        }
    }
  return aBeginVoxel;
}




