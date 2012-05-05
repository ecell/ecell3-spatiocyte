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
  for(int i(0); i != 9; ++i)
    {
      Voxel* aVoxel(getBeginVoxel());
      theProcessSpecies[i]->addMolecule(aVoxel);
      Point S(theSpatiocyteStepper->coord2point(aVoxel->coord));
      std::cout << "S.x:" << S.x << " y:" << S.y << " z:" << S.z << std::endl;
      Voxel* aNeighbor(NULL);
      do {
        aNeighbor = getNeighbor(aVoxel, S); 
        if(aNeighbor)
          {
            std::cout << "S.x:" << S.x << " y:" << S.y << " z:" << S.z << std::endl;
            theProcessSpecies[i]->addMolecule(aNeighbor);
            aVoxel = aNeighbor;
          }
      } while(aNeighbor);
    }
  theProcessSpecies[0]->setIsPopulated();
}

Voxel* OneDCompartmentProcess::getNeighbor(Voxel* aVoxel, Point& S)
{
  double shortestDist(1e+10);
  Voxel* aNeighbor(NULL);
  Point currS(S);
  for(int j(0); j != aVoxel->adjoiningSize; ++j)
    {
      Point N(theSpatiocyteStepper->coord2point(
                                aVoxel->adjoiningVoxels[j]->coord));
      if(aVoxel->adjoiningVoxels[j]->id == theComp->vacantID)
        {
          //std::cout << "N.x:" << N.x << " y:" << N.y << " z:" << N.z << std::endl;
          double t((-E.x*N.x-E.y*N.y-E.z*N.z+E.x*S.x+E.y*S.y+E.z*S.z+N.x*W.x-S.x*W.x+N.y*W.y-S.y*W.y+N.z*W.z-S.z*W.z)/(E.x*E.x+E.y*E.y+E.z*E.z-2*E.x*W.x+W.x*W.x-2*E.y*W.y+W.y*W.y-2*E.z*W.z+W.z*W.z));
          //std::cout << "t:" << t << std::endl;
          if(t<0)
            {
              double dist(sqrt(pow(-N.x+S.x+t*(-E.x+W.x),2)+pow(-N.y+S.y+t*(-E.y+W.y),2)+pow(-N.z+S.z+t*(-E.z+W.z),2)));
              if(dist < shortestDist && dist < 0.75)
                {
                  std::cout << "dist:" << dist << std::endl;
                  aNeighbor = aVoxel->adjoiningVoxels[j];
                  shortestDist = dist;
                  currS.x = S.x+t*(-E.x+W.x);
                  currS.y = S.y+t*(-E.y+W.y);
                  currS.z = S.z+t*(-E.z+W.z);
                }
            }
        }
    } 
  S = currS;
  return aNeighbor;
}

Voxel* OneDCompartmentProcess::getBeginVoxel()
{
  Point C(theComp->centerPoint);
  //East point
  E.x = theComp->lengthX/2;
  E.y = 0;
  E.z = 0;
  //std::cout << "E.x:" << E.x << " y:" << E.y << " z:" << E.z << std::endl;
  //West point
  W.x = -theComp->lengthX/2;
  W.y = 0;
  W.z = 0;
  //std::cout << "W.x:" << W.x << " y:" << W.y << " z:" << W.z << std::endl;
  theSpatiocyteStepper->rotateX(theComp->rotateX, &E);
  theSpatiocyteStepper->rotateY(theComp->rotateY, &E);
  theSpatiocyteStepper->rotateZ(theComp->rotateZ, &E);
  theSpatiocyteStepper->rotateX(theComp->rotateX, &W);
  theSpatiocyteStepper->rotateY(theComp->rotateY, &W);
  theSpatiocyteStepper->rotateZ(theComp->rotateZ, &W);
  //std::cout << "E.x:" << E.x << " y:" << E.y << " z:" << E.z << std::endl;
  //std::cout << "W.x:" << W.x << " y:" << W.y << " z:" << W.z << std::endl;
  //std::cout << "C.x:" << C.x << " y:" << C.y << " z:" << C.z << std::endl;
  E.x += C.x;
  E.y += C.y;
  E.z += C.z;
  W.x += C.x;
  W.y += C.y;
  W.z += C.z;
  std::cout << "E.x:" << E.x << " y:" << E.y << " z:" << E.z << std::endl;
  std::cout << "W.x:" << W.x << " y:" << W.y << " z:" << W.z << std::endl;
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




