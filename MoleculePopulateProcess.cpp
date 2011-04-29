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
      aSpecies->setIsPopulated();
    }
}

void MoleculePopulateProcess::populateUniformSparse(Species* aSpecies)
{
  Comp* aComp(aSpecies->getComp());
  if(!aSpecies->getIsPopulated())
    {
      unsigned int aSize(aSpecies->getPopulateMoleculeSize());
      if(UniformRadius == 1 && !OriX && !OriY && !OriZ)
        {
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
          double diaX(UniformRadius*aComp->lengthX*(1-sqrt(OriX*OriX)));
          double diaY(UniformRadius*aComp->lengthY*(1-sqrt(OriY*OriY)));
          double diaZ(UniformRadius*aComp->lengthZ*(1-sqrt(OriZ*OriZ)));
          double dia(diaX);
          if(diaY < dia)
            {
              dia = diaY;
            }
          if(diaZ < dia)
            {
              dia = diaZ;
            }
          Point aPoint;
          aPoint.x = aComp->lengthX*OriX*0.5+aComp->centerPoint.x;
          aPoint.y = aComp->lengthY*OriY*0.5+aComp->centerPoint.y;
          aPoint.z = aComp->lengthZ*OriZ*0.5+aComp->centerPoint.z;
          std::cout << "the Comp center point:" << std::endl;
          std::cout << "x:" << aComp->centerPoint.x <<
            " y:" << aComp->centerPoint.y <<
            " z:" << aComp->centerPoint.z << std::endl;
          std::cout << "the placement target point:" << std::endl;
          std::cout << "x:" << aPoint.x << " y:" << aPoint.y << " z:" << aPoint.z
            << std::endl;
          for(unsigned int j(0); j != aSize; ++j)
            {
              Voxel* aVoxel;
              do
                {
                  double aRandX(gsl_rng_uniform(getStepper()->getRng()));
                  double aRandY(gsl_rng_uniform(getStepper()->getRng()));
                  double aRandZ(gsl_rng_uniform(getStepper()->getRng()));
                  Point aTargetPoint;
                  aTargetPoint.x = aPoint.x+dia*(aRandX-0.5); 
                  aTargetPoint.y = aPoint.y+dia*(aRandY-0.5); 
                  aTargetPoint.z = aPoint.z+dia*(aRandZ-0.5);
                  aVoxel = theSpatiocyteStepper->point2voxel(aTargetPoint);
                }
              while(aVoxel->id != aComp->vacantID);
              aSpecies->addMolecule(aVoxel);
            }
        }
      aSpecies->setIsPopulated();
    }
}

