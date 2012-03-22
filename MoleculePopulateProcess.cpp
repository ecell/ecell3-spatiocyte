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

#include <algorithm>
#include <gsl/gsl_randist.h>
#include <boost/lexical_cast.hpp>
#include "MoleculePopulateProcess.hpp"
#include "SpatiocyteSpecies.hpp"

LIBECS_DM_INIT(MoleculePopulateProcess, Process); 

void MoleculePopulateProcess::initializeSecond()
{
  SpatiocyteProcess::initializeSecond();
  for(std::vector<Species*>::const_iterator i(theProcessSpecies.begin());
      i != theProcessSpecies.end(); ++i)
    {
      (*i)->setPopulateProcess(this, GaussianSigma);
    }
}

void MoleculePopulateProcess::fire()
{
  for(std::vector<Species*>::const_iterator i(theProcessSpecies.begin());
      i != theProcessSpecies.end(); ++i)
    {
      (*i)->removeMolecules();
      populateUniformSparse(*i);
    }
  theStepInterval = ResetTime;
  theTime += theStepInterval; 
  thePriorityQueue->move(theQueueID);
}

void MoleculePopulateProcess::populateGaussian(Species* aSpecies)
{
}

void MoleculePopulateProcess::populateUniformDiffuseVacant(Species* aSpecies)
{
  std::cout << "   Populating:" << getIDString(aSpecies) << std::endl;
  if(!aSpecies->getIsPopulated())
    {
      if(UniformRadiusX == 1 && UniformRadiusY == 1 && UniformRadiusZ == 1 &&
         !OriginX && !OriginY && !OriginZ)
        {
          Species* aVacantSpecies(aSpecies->getVacantSpecies());
          unsigned int aSize(aSpecies->getPopulateMoleculeSize());
          for(unsigned int i(0); i != aSize; ++i)
            {
              Voxel* aMolecule(aVacantSpecies->getRandomMolecule());
              aSpecies->addMolecule(aMolecule);
            }
        }
      else
        {
          populateUniformRanged(aSpecies);
        }
      aSpecies->setIsPopulated();
    }
}

void MoleculePopulateProcess::populateUniformDense(Species* aSpecies,
                                              unsigned int aList[], 
                                              unsigned int* aCount)
{
  std::cout << "   Populating:" << getIDString(aSpecies) << std::endl;
  Comp* aComp(aSpecies->getComp());
  if(!aSpecies->getIsPopulated())
    {
      if(UniformRadiusX == 1 && UniformRadiusY == 1 && UniformRadiusZ == 1 &&
         !OriginX && !OriginY && !OriginZ)
        {
          unsigned int aSize(aSpecies->getPopulateMoleculeSize());
          for(unsigned int j(0); j != aSize; ++j)
            {
              Voxel* aVoxel;
              do
                {
                  aVoxel = theSpatiocyteStepper->coord2voxel(
                    aComp->coords[aList[(*aCount)++]]);
                }
              while(aVoxel->id != aComp->vacantID);
              aSpecies->addMolecule(aVoxel);
            }
        }
      else
        {
          populateUniformRanged(aSpecies);
        }
      aSpecies->setIsPopulated();
    }
}

void MoleculePopulateProcess::populateUniformSparse(Species* aSpecies)
{
  Comp* aComp(aSpecies->getComp());
  if(!aSpecies->getIsPopulated())
    {
      if(UniformRadiusX == 1 && UniformRadiusY == 1 && UniformRadiusZ == 1 &&
         !OriginX && !OriginY && !OriginZ)
        {
          unsigned int aSize(aSpecies->getPopulateMoleculeSize());
          int availableVoxelSize(aComp->coords.size());
          for(unsigned int j(0); j != aSize; ++j)
            {
              Voxel* aVoxel;
              do
                {
                  aVoxel = theSpatiocyteStepper->coord2voxel(
                     aComp->coords[gsl_rng_uniform_int(
                                getStepper()->getRng(), availableVoxelSize)]);
                }
              while(aVoxel->id != aComp->vacantID);
              aSpecies->addMolecule(aVoxel);
            }
        }
      else
        { 
          populateUniformRanged(aSpecies);
        }
      aSpecies->setIsPopulated();
    }
}

void MoleculePopulateProcess::populateUniformRanged(Species* aSpecies)
{
  Comp* aComp(aSpecies->getComp());
  double delta(0);
  // Increase the compartment dimensions by delta if it is a surface 
  // compartment:
  if(aComp->dimension == 2)
    {
      delta = 0.1;
    }
  double maxX(std::min(1.0, OriginX+UniformRadiusX));
  double minX(std::max(-1.0, OriginX-UniformRadiusX));
  double maxY(std::min(1.0, OriginY+UniformRadiusY));
  double minY(std::max(-1.0, OriginY-UniformRadiusY));
  double maxZ(std::min(1.0, OriginZ+UniformRadiusZ));
  double minZ(std::max(-1.0, OriginZ-UniformRadiusZ)); 
  maxX = aComp->centerPoint.x + maxX*aComp->lengthX/2*(1+delta);
  minX = aComp->centerPoint.x + minX*aComp->lengthX/2*(1+delta);
  maxY = aComp->centerPoint.y + maxY*aComp->lengthY/2*(1+delta);
  minY = aComp->centerPoint.y + minY*aComp->lengthY/2*(1+delta);
  maxZ = aComp->centerPoint.z + maxZ*aComp->lengthZ/2*(1+delta);
  minZ = aComp->centerPoint.z + minZ*aComp->lengthY/2*(1+delta);
  std::vector<unsigned int> aCoords;
  for(std::vector<unsigned int>::iterator i(aComp->coords.begin());
      i != aComp->coords.end(); ++i)
    {
      Voxel* aVoxel(theSpatiocyteStepper->coord2voxel(*i));
      Point aPoint(theSpatiocyteStepper->coord2point(aVoxel->coord));
      if(aVoxel->id == aSpecies->getVacantID() &&
         aPoint.x < maxX && aPoint.x > minX &&
         aPoint.y < maxY && aPoint.y > minY &&
         aPoint.z < maxZ && aPoint.z > minZ)
        {
          aCoords.push_back(*i);
        }
    }
  unsigned int aSize(aSpecies->getPopulateMoleculeSize());
  if(aCoords.size() < aSize)
    {
      THROW_EXCEPTION(ValueError, String(
                      getPropertyInterface().getClassName()) +
                      "[" + getFullID().asString() + "]: The are " +
                      int2str(aSize) + " " + getIDString(aSpecies) +
                      " molecules that must be populated,\n but there are " +
                      "only " + int2str(aCoords.size()) + " vacant voxels in" +
                      getIDString(aComp) + " that can be populated.");
    }
  unsigned int aCoordsArray[aCoords.size()]; 
  for(unsigned int i(0); i != aCoords.size(); ++i)
    {
      aCoordsArray[i] = aCoords[i];
    }
  gsl_ran_shuffle(getStepper()->getRng(), aCoordsArray, aCoords.size(),
                  sizeof(unsigned int));
  for(unsigned int i(0); i != aSize; ++i)
    {
      aSpecies->addMolecule(theSpatiocyteStepper->coord2voxel(aCoordsArray[i]));
    }
}
