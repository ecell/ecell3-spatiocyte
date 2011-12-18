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

#include "MoleculePopulateProcess.hpp"
#include "SpatiocyteSpecies.hpp"
#include <algorithm>

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

void MoleculePopulateProcess::populateUniformDense(Species* aSpecies,
                                              unsigned int aList[], 
                                              unsigned int* aCount)
{
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
  unsigned int aSize(aSpecies->getPopulateMoleculeSize());
  for(unsigned int j(0); j != aSize; ++j)
    {
      Voxel* aVoxel;
      do
        {
          double aRandX(gsl_rng_uniform(getStepper()->getRng()));
          double aRandY(gsl_rng_uniform(getStepper()->getRng()));
          double aRandZ(gsl_rng_uniform(getStepper()->getRng()));
          Point aTargetPoint;
          aTargetPoint.x = aComp->centerPoint.x +
            (minX+(maxX-minX)*aRandX)*aComp->lengthX/2*(1+delta);
          aTargetPoint.y = aComp->centerPoint.y + 
            (minY+(maxY-minY)*aRandY)*aComp->lengthY/2*(1+delta);
          aTargetPoint.z = aComp->centerPoint.z + 
            (minZ+(maxZ-minZ)*aRandZ)*aComp->lengthZ/2*(1+delta);
          aVoxel = theSpatiocyteStepper->point2voxel(aTargetPoint);
        }
      while(aVoxel->id != aComp->vacantID);
      aSpecies->addMolecule(aVoxel);
    }
}

