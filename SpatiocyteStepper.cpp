//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of E-Cell Simulation Environment package
//
//                Copyright (C) 2006-2009 Keio University
//
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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


#include <time.h>
#include <gsl/gsl_randist.h>
#include <libecs/Model.hpp>
#include <libecs/System.hpp>
#include <libecs/Process.hpp>
#include "SpatiocyteStepper.hpp"
#include "SpatiocyteSpecies.hpp"
#include "SpatiocyteProcessInterface.hpp"
#include "ReactionProcessInterface.hpp"

LIBECS_DM_INIT(SpatiocyteStepper, Stepper);

void SpatiocyteStepper::initialize()
{
  if(isInitialized)
    {
      return;
    }
  isInitialized = true;
  Stepper::initialize(); 
  if(theProcessVector.empty())
    {
      THROW_EXCEPTION(InitializationFailed,
                      getPropertyInterface().getClassName() + 
                      ": at least one Process must be defined in this" +
                      " Stepper.");
    } 
  //We need a Comp tree to assign the voxels to each Comp
  //and get the available number of vacant voxels. The compartmentalized
  //vacant voxels are needed to randomly place molecules according to the
  //Comp:
  std::cout << "1. creating compartments..." << std::endl;
  registerComps();
  setCompsProperties();
  std::cout << "2. setting up lattice properties..." << std::endl;
  setLatticeProperties(); 
  setCompsCenterPoint();
  //All species have been created at this point, we initialize them now:
  std::cout << "3. initializing species..." << std::endl;
  initSpecies();
  std::cout << "4. initializing processes the second time..." << std::endl;
  initProcessSecond();
  std::cout << "5. constructing lattice..." << std::endl;
  constructLattice();
  std::cout << "6. setting intersecting compartment list..." << std::endl;
  setIntersectingCompartmentList();
  std::cout << "7. compartmentalizing lattice..." << std::endl;
  compartmentalizeLattice();
  std::cout << "8. setting up compartment voxels properties..." << std::endl;
  setCompVoxelProperties();
  std::cout << "9. printing simulation parameters..." << std::endl;
  storeSimulationParameters();
  printSimulationParameters();
  std::cout << "10. populating compartments with molecules..." << std::endl;
  populateComps();
  std::cout << "11. initializing processes the third time..." << std::endl;
  initProcessThird();
  std::cout << "12. initializing the priority queue..." << std::endl;
  initPriorityQueue();
  std::cout << "13. initializing processes the fourth time..." << std::endl;
  initProcessFourth();
  std::cout << "14. initializing processes the last time..." << std::endl;
  initProcessLastOnce();
  std::cout << "15. printing final process parameters..." << std::endl;
  printProcessParameters();
  std::cout << "16. simulation is started..." << std::endl;
}

unsigned int SpatiocyteStepper::getStartCoord()
{
  return theStartCoord;
}

unsigned int SpatiocyteStepper::getRowSize()
{
  return theRowSize;
}

unsigned int SpatiocyteStepper::getLayerSize()
{
  return theLayerSize;
}

unsigned int SpatiocyteStepper::getColSize()
{
  return theColSize;
}

unsigned int SpatiocyteStepper::getLatticeSize()
{
  return theLattice.size();
}

Point SpatiocyteStepper::getCenterPoint()
{
  return theCenterPoint;
} 

double SpatiocyteStepper::getNormalizedVoxelRadius()
{
  return theNormalizedVoxelRadius;
}

void SpatiocyteStepper::reset(int seed)
{
  gsl_rng_set(getRng(), seed); 
  setCurrentTime(0);
  initProcessSecond();
  clearComps();
  populateComps();
  initProcessThird();
  initPriorityQueue();
  initProcessFourth();
  //checkLattice();
}

Species* SpatiocyteStepper::addSpecies(Variable* aVariable)
{
  std::vector<Species*>::iterator aSpeciesIter(variable2ispecies(aVariable));
  if(aSpeciesIter == theSpecies.end())
    {
      Species *aSpecies(new Species(this, aVariable, theSpecies.size(),
                                    (int)aVariable->getValue(), getRng()));
      theSpecies.push_back(aSpecies);
      return aSpecies;
    }
  return *aSpeciesIter;
}

Species* SpatiocyteStepper::getSpecies(Variable* aVariable)
{
  std::vector<Species*>::iterator aSpeciesIter(variable2ispecies(aVariable));
  if(aSpeciesIter == theSpecies.end())
    {
      return NULL;
    }
  return *aSpeciesIter;
}

std::vector<Species*> SpatiocyteStepper::getSpecies()
{
  return theSpecies;
}

bool SpatiocyteStepper::isBoundaryCoord(unsigned int aCoord, bool isVolume)
{
  //This method is only for checking boundaries on a cube:
  unsigned int aRow;
  unsigned int aLayer;
  unsigned int aCol;
  coord2global(aCoord, &aRow, &aLayer, &aCol);
  if(isVolume)
    {
      //If the voxel is on one of the 6 cuboid surfaces:
      if(aRow == 1 || aLayer == 1 || aCol == 1 ||
         aRow == theRowSize-2 || aLayer == theLayerSize-2 || 
         aCol == theColSize-2)
        {
          return true;
        }
    }
  //If it is a surface voxel:
  else
    {
      //If the voxel is on one of the 12 edges of the cube:
      if((aRow <= 1 && aCol <= 1) ||
         (aRow <= 1 && aLayer <= 1) ||
         (aCol <= 1 && aLayer <= 1) ||
         (aRow <= 1 && aCol >= theColSize-2) ||
         (aRow <= 1 && aLayer >= theLayerSize-2) ||
         (aCol <= 1 && aRow >= theRowSize-2) ||
         (aCol <= 1 && aLayer >= theLayerSize-2) ||
         (aLayer <= 1 && aRow >= theRowSize-2) ||
         (aLayer <= 1 && aCol >= theColSize-2) ||
         (aRow >= theRowSize-2 && aCol >= theColSize-2) ||
         (aRow >= theRowSize-2 && aLayer >= theLayerSize-2) ||
         (aCol >= theColSize-2 && aLayer >= theLayerSize-2))
        {
          return true;
        }
    }
  return false;
}

Voxel* SpatiocyteStepper::getPeriodicVoxel(unsigned int aCoord,
                                           bool isVolume,
                                           Origin* anOrigin)
{
  //This method is only for checking boundaries on a cube:
  unsigned int aRow;
  unsigned int aLayer;
  unsigned int aCol;
  coord2global(aCoord, &aRow, &aLayer, &aCol);
  unsigned int nextRow(aRow);
  unsigned int nextLayer(aLayer);
  unsigned int nextCol(aCol);
  unsigned int adj(1);
  if(isVolume)
    {
      adj = 0;
    }
  if(aRow == 1+adj)
    {
      nextRow = theRowSize-(3+adj);
      --anOrigin->row;
    }
  else if(aRow == theRowSize-(2+adj))
    {
      nextRow = 2+adj;
      ++anOrigin->row;
    }
  if(aLayer == 1+adj)
    {
      nextLayer = theLayerSize-(3+adj);
      --anOrigin->layer;
    }
  else if(aLayer == theLayerSize-(2+adj))
    {
      nextLayer = 2+adj;
      ++anOrigin->layer;
    }
  if(aCol == 1+adj)
    {
      nextCol = theColSize-(3+adj);
      --anOrigin->col;
    }
  else if(aCol == theColSize-(2+adj))
    {
      nextCol = 2+adj;
      ++anOrigin->col;
    }
  if(nextRow != aRow || nextCol != aCol || nextLayer != aLayer)
    {
      return &theLattice[nextRow+
                         theRowSize*nextLayer+
                         theRowSize*theLayerSize*nextCol];
    }
  return NULL;
}

Point SpatiocyteStepper::getPeriodicPoint(unsigned int aCoord,
                                          bool isVolume,
                                          Origin* anOrigin)
{
  unsigned int adj(1);
  if(isVolume)
    {
      adj = 0;
    }
  unsigned int row(theRowSize-(3+adj)-(1+adj));
  unsigned int layer(theLayerSize-(3+adj)-(1+adj));
  unsigned int col(theColSize-(3+adj)-(1+adj));
  unsigned int aGlobalCol;
  unsigned int aGlobalLayer;
  unsigned int aGlobalRow;
  coord2global(aCoord, &aGlobalRow, &aGlobalLayer, &aGlobalCol);
  int aRow(aGlobalRow+anOrigin->row*row);
  int aLayer(aGlobalLayer+anOrigin->layer*layer);
  int aCol(aGlobalCol+anOrigin->col*col);
  Point aPoint;
  switch(LatticeType)
    {
    case HCP_LATTICE: 
      aPoint.y = (aCol%2)*theHCPk+theHCPl*aLayer;
      aPoint.z = aRow*2*theNormalizedVoxelRadius+
        ((aLayer+aCol)%2)*theNormalizedVoxelRadius;
      aPoint.x = aCol*theHCPh;
      break;
    case CUBIC_LATTICE:
      aPoint.y = aLayer*2*theNormalizedVoxelRadius;
      aPoint.z = aRow*2*theNormalizedVoxelRadius;
      aPoint.x = aCol*2*theNormalizedVoxelRadius;
      break;
    }
  return aPoint;
}


std::vector<Species*>::iterator
SpatiocyteStepper::variable2ispecies(Variable* aVariable)
{
  for(std::vector<Species*>::iterator i(theSpecies.begin());
      i != theSpecies.end(); ++i)
    {
      if((*i)->getVariable() == aVariable)
        {
          return i;
        }
    }
  return theSpecies.end();
} 

Species* SpatiocyteStepper::variable2species(Variable* aVariable)
{
  for(std::vector<Species*>::iterator i(theSpecies.begin());
      i != theSpecies.end(); ++i)
    {
      if((*i)->getVariable() == aVariable)
        {
          return (*i);
        }
    }
  return NULL;
}

void SpatiocyteStepper::checkLattice()
{
  std::vector<int> list;
  for(unsigned int i(0); i!=theSpecies.size(); ++i)
    {
      list.push_back(0);
    }
  for(unsigned int i(0); i!=theLattice.size(); ++i)
    {
      ++list[theLattice[i].id];
    }
  int volumeCnt(0);
  int surfaceCnt(0);
  for(unsigned int i(0); i!=list.size(); ++i)
    {
      std::cout << "i:" << i << " ";
      if(theSpecies[i]->getVariable() != NULL)
        {
          std::cout << theSpecies[i]->getVariable()->getFullID().asString();
          if(theSpecies[i]->getIsVolume())
            {
              volumeCnt += list[i];
            }
          else
            {
              surfaceCnt += list[i];
            }
        }
      std::cout << " cnt:" << list[i] << std::endl;
    }
  std::cout << "total volume:" << volumeCnt << std::endl;
  std::cout << "total surface:" << surfaceCnt << std::endl;
  std::cout << "total volume+surface:" << surfaceCnt+volumeCnt << std::endl;
}

void SpatiocyteStepper::checkSurfaceComp()
{
  for(unsigned int i(0); i!=theLattice.size(); ++i)
    {
      Voxel* aVoxel(&theLattice[i]);
      for(unsigned int k(0); k != theAdjoiningVoxelSize; ++k)
        {
          if(aVoxel == aVoxel->adjoiningVoxels[k])
            {
              theSpecies[1]->addSimpleMolecule(aVoxel);
              break;
            }
        }
    }
}

void SpatiocyteStepper::initSpecies()
{
  for(std::vector<Species*>::iterator i(theSpecies.begin());
      i != theSpecies.end(); ++i)
    {
      (*i)->initialize(theSpecies.size(), theAdjoiningVoxelSize);
    }
}

void SpatiocyteStepper::initProcessSecond()
{
  for(std::vector<Process*>::const_iterator i(theProcessVector.begin());
      i != theProcessVector.end(); ++i)
    {      
      SpatiocyteProcessInterface*
        aProcess(dynamic_cast<SpatiocyteProcessInterface*>(*i));
      aProcess->initializeSecond();
    }
}

void SpatiocyteStepper::printProcessParameters()
{
  for(std::vector<Process*>::const_iterator i(theProcessVector.begin());
      i != theProcessVector.end(); ++i)
    {      
      SpatiocyteProcessInterface*
        aProcess(dynamic_cast<SpatiocyteProcessInterface*>(*i));
      aProcess->printParameters();
    }
  std::cout << std::endl;
}

void SpatiocyteStepper::initProcessThird()
{
  for(std::vector<Process*>::const_iterator i(theProcessVector.begin());
      i != theProcessVector.end(); ++i)
    {      
      SpatiocyteProcessInterface*
        aProcess(dynamic_cast<SpatiocyteProcessInterface*>(*i));
      aProcess->initializeThird();
    }
}

void SpatiocyteStepper::initProcessFourth()
{
  for(std::vector<Process*>::const_iterator i(theProcessVector.begin());
      i != theProcessVector.end(); ++i)
    {      
      SpatiocyteProcessInterface*
        aProcess(dynamic_cast<SpatiocyteProcessInterface*>(*i));
      aProcess->initializeFourth();
    }
  setStepInterval(thePriorityQueue.getTop()->getTime()-getCurrentTime());
}

void SpatiocyteStepper::initProcessLastOnce()
{
  for(std::vector<Process*>::const_iterator i(theProcessVector.begin());
      i != theProcessVector.end(); ++i)
    {      
      SpatiocyteProcessInterface*
        aProcess(dynamic_cast<SpatiocyteProcessInterface*>(*i));
      aProcess->initializeLastOnce();
    }
}

void SpatiocyteStepper::initPriorityQueue()
{
  const double aCurrentTime(getCurrentTime());
  thePriorityQueue.clear();
  for(std::vector<Process*>::const_iterator i(theProcessVector.begin());
      i != theProcessVector.end(); ++i)
    {      
      Process* const aProcess(*i);
      SpatiocyteProcessInterface*
        aSpatiocyteProcess(dynamic_cast<SpatiocyteProcessInterface*>(*i));
      if(aSpatiocyteProcess != NULL)
        {
          aSpatiocyteProcess->setTime(aCurrentTime+aProcess->getStepInterval());
          aSpatiocyteProcess->setPriorityQueue(&thePriorityQueue);
          String aClassName(aProcess->getPropertyInterface().getClassName());
          //Not all SpatiocyteProcesses are inserted into the priority queue.
          //Only the following processes are inserted in the PriorityQueue and
          //executed at simulation steps according to their execution times:
          if(aClassName == "DiffusionProcess" ||
             aClassName == "IteratingLogProcess" ||
             aClassName == "CompartmentGrowthProcess" ||
             aClassName == "MoleculePopulateProcess" ||
             aClassName == "CoordinateLogProcess" ||
             aClassName == "VisualizationLogProcess" ||
             aClassName == "H5VisualizationLogProcess" ||
             aClassName == "MicroscopyTrackingProcess" ||
             aClassName == "OscillationAnalysisProcess" ||
             aClassName == "PeriodicBoundaryDiffusionProcess" ||
             aClassName == "SpatiocyteNextReactionProcess" ||
             aClassName == "PolymerFragmentationProcess")
            {
              aSpatiocyteProcess->setQueueID(
                                   thePriorityQueue.push(aSpatiocyteProcess));
            }
        }
      //The following processes never interrupt other Processes.
      //We exclude them here and set up the interrupt for the remaining
      //processes. All processes which interrupt other processes have
      //the ReactionProcess as the base class.
      ReactionProcessInterface*
        aReactionProcess(dynamic_cast<ReactionProcessInterface*>(*i));
      if(aReactionProcess != NULL)
        {
          aReactionProcess->setInterrupt(theProcessVector, *i);
        }
    } 
}

void SpatiocyteStepper::populateComps()
{
  for(std::vector<Comp*>::const_iterator i(theComps.begin());
      i != theComps.end(); ++i)
    {
      populateComp(*i);
    }
}

void SpatiocyteStepper::clearComps()
{
  for(std::vector<Comp*>::const_iterator i(theComps.begin());
      i != theComps.end(); ++i)
    {
      clearComp(*i);
    }
}


inline void SpatiocyteStepper::step()
{
  do
    {
      thePriorityQueue.getTop()->fire();
    }
  while(thePriorityQueue.getTop()->getTime() == getCurrentTime());
  setNextTime(thePriorityQueue.getTop()->getTime());
  //checkLattice();
} 


void SpatiocyteStepper::registerComps()
{
  System* aRootSystem(getModel()->getRootSystem());
  std::vector<Comp*> allSubs;
  //The root Comp is theComps[0]
  theComps.push_back(registerComp(aRootSystem, &allSubs));
  //After this we will create an species to get an ID to represent
  //NULL Comps. So let us specify the correct
  //size of the biochemical species to be simulated before additional
  //non-biochemical species are created:
  theBioSpeciesSize = theSpecies.size();
  //Create one last species to represent a NULL Comp. This is for
  //voxels that do not belong to any Comps:
  Species* aSpecies(new Species(this, NULL, theSpecies.size(), 0, getRng()));
  theSpecies.push_back(aSpecies);
  aSpecies->setComp(NULL);
  theNullID = aSpecies->getID(); 
  //Expand the tree of immediate subComps into single list such that
  //the super Comps come first while the subComps 
  //come later in the list:
  std::vector<Comp*> Comps(theComps[0]->immediateSubs);
  while(!Comps.empty())
    {
      std::vector<Comp*> subComps;
      for(unsigned int i(0); i != Comps.size(); ++i)
        {
          theComps.push_back(Comps[i]);
          for(unsigned int j(0);
              j != Comps[i]->immediateSubs.size(); ++j)
            {
              subComps.push_back(Comps[i]->immediateSubs[j]);
            }
        }
      Comps = subComps;
    }
}

//allSubs contains all the subComps (child, grand child, great grand
//child, etc). Used to calculate the total number of Comp voxels.
Comp* SpatiocyteStepper::registerComp(System* aSystem,
                                            std::vector<Comp*>* allSubs)
{ 
  //We execute this function to register the System, and its subsystems
  //recursively.
  Comp* aComp(new Comp);
  aComp->minRow = UINT_MAX;
  aComp->minCol = UINT_MAX;
  aComp->minLayer = UINT_MAX;
  aComp->maxRow = 0;
  aComp->maxCol = 0;
  aComp->maxLayer = 0;
  aComp->lengthX = 0;
  aComp->lengthY = 0;
  aComp->lengthZ = 0;
  aComp->originX = 0;
  aComp->originY = 0;
  aComp->originZ = 0;
  aComp->rotateX = 0;
  aComp->rotateY = 0;
  aComp->rotateZ = 0;
  aComp->xyPlane = 0;
  aComp->xzPlane = 0;
  aComp->yzPlane = 0;
  aComp->specVolume = 0;
  aComp->system = aSystem;
  aComp->surfaceSub = NULL;
  aComp->diffusiveComp = NULL;
  //Default Comp geometry is Cuboid:
  aComp->geometry = 0;
  //Default is volume Comp:
  aComp->dimension = 3;
  if(getVariable(aSystem, "DIMENSION"))
    { 
      const double aDimension(aSystem->getVariable("DIMENSION")->getValue());
      aComp->dimension = aDimension;
    }
  if(aComp->dimension == 3)
    {
      if(getVariable(aSystem, "GEOMETRY"))
        { 
          aComp->geometry = aSystem->getVariable("GEOMETRY")->getValue();
        }
      if(getVariable(aSystem, "LENGTHX"))
        {
          aComp->lengthX = aSystem->getVariable("LENGTHX")->getValue();
        }
      if(getVariable(aSystem, "LENGTHY"))
        {
          aComp->lengthY = aSystem->getVariable("LENGTHY")->getValue();
        }
      if(getVariable(aSystem, "LENGTHZ"))
        {
          aComp->lengthZ = aSystem->getVariable("LENGTHZ")->getValue();
        }
      if(getVariable(aSystem, "ORIGINX"))
        {
          aComp->originX = aSystem->getVariable("ORIGINX")->getValue();
        }
      if(getVariable(aSystem, "ORIGINY"))
        {
          aComp->originY = aSystem->getVariable("ORIGINY")->getValue();
        }
      if(getVariable(aSystem, "ORIGINZ"))
        {
          aComp->originZ = aSystem->getVariable("ORIGINZ")->getValue();
        }
      if(getVariable(aSystem, "ROTATEX"))
        {
          aComp->rotateX = aSystem->getVariable("ROTATEX")->getValue();
        }
      if(getVariable(aSystem, "ROTATEY"))
        {
          aComp->rotateY = aSystem->getVariable("ROTATEY")->getValue();
        }
      if(getVariable(aSystem, "ROTATEZ"))
        {
          aComp->rotateZ = aSystem->getVariable("ROTATEZ")->getValue();
        }
      if(getVariable(aSystem, "XYPLANE"))
        {
          aComp->xyPlane = aSystem->getVariable("XYPLANE")->getValue();
        }
      if(getVariable(aSystem, "XZPLANE"))
        {
          aComp->xzPlane = aSystem->getVariable("XZPLANE")->getValue();
        }
      if(getVariable(aSystem, "YZPLANE"))
        {
          aComp->yzPlane = aSystem->getVariable("YZPLANE")->getValue();
        }
    }
  registerCompSpecies(aComp);
  //Systems contains all the subsystems of a System.
  //For example /membrane is the subsystem of /:
  FOR_ALL(System::Systems, aSystem->getSystems())
    {
      std::vector<Comp*> currSubs; 
      Comp* aSubComp(registerComp(i->second, &currSubs)); 
      allSubs->push_back(aSubComp);
      for(unsigned int j(0); j != currSubs.size(); ++j)
        {
          allSubs->push_back(currSubs[j]);
        }
      aComp->immediateSubs.push_back(aSubComp);
      if(aSubComp->dimension == 2)
      {
          aSubComp->geometry = aComp->geometry;
          aSubComp->specVolume = aComp->specVolume;
          aSubComp->lengthX = aComp->lengthX;
          aSubComp->lengthY = aComp->lengthY;
          aSubComp->lengthZ = aComp->lengthZ;
          aSubComp->originX = aComp->originX;
          aSubComp->originY = aComp->originY;
          aSubComp->originZ = aComp->originZ;
          aSubComp->xyPlane = aComp->xyPlane;
          aSubComp->xzPlane = aComp->xzPlane;
          aSubComp->yzPlane = aComp->yzPlane;

          for(unsigned int i(0); i != aSubComp->lineSubs.size(); ++i)
          {
              Comp* lineComp(aSubComp->lineSubs[i]);
              lineComp->geometry = aComp->geometry;
              lineComp->specVolume = aComp->specVolume;
              lineComp->lengthX = aComp->lengthX;
              lineComp->lengthY = aComp->lengthY;
              lineComp->lengthZ = aComp->lengthZ;
              lineComp->originX = aComp->originX;
              lineComp->originY = aComp->originY;
              lineComp->originZ = aComp->originZ;
              lineComp->xyPlane = aComp->xyPlane;
              lineComp->xzPlane = aComp->xzPlane;
              lineComp->yzPlane = aComp->yzPlane;
          }

          aComp->surfaceSub = aSubComp;
      }
      else if(aSubComp->dimension == 1)
      {
          // Properties of aComp (dimension == 2) is not fully set here yet.
          aComp->lineSubs.push_back(aSubComp);
      }
    }
  aComp->allSubs = *allSubs;
  return aComp;
}

void SpatiocyteStepper::setCompsProperties()
{
  for(unsigned int i(0); i != theComps.size(); ++i)
    {
      setCompProperties(theComps[i]);
    }
}

void SpatiocyteStepper::setCompsCenterPoint()
{
  for(unsigned int i(0); i != theComps.size(); ++i)
    {
      setCompCenterPoint(theComps[i]);
    }
}

void SpatiocyteStepper::registerCompSpecies(Comp* aComp)
{
  System* aSystem(aComp->system);
  FOR_ALL(System::Variables, aSystem->getVariables())
    {
      Variable* aVariable(i->second);
      if(aVariable->getID() == "VACANT")
        {
          aComp->enclosed = aVariable->getValue();
          //Set the number of vacant molecules to be always 0 because
          //when we populate lattice we shouldn't create more vacant
          //molecules than the ones already created for the Comp:
          aVariable->setValue(0);
          Species* aSpecies(addSpecies(aVariable));
          aComp->vacantID = aSpecies->getID();
          aSpecies->setIsVacant();
        }
      std::vector<Species*>::iterator j(variable2ispecies(aVariable));
      if(j != theSpecies.end())
        {
          aComp->species.push_back(*j);
          (*j)->setComp(aComp);
          if(aComp->dimension == 3)
            {
              //By default all biochemical species are assumed to be
              //surface species. Here we set it as a volume species:
              (*j)->setIsVolume();
            }
        }
    }
}

void SpatiocyteStepper::setLatticeProperties()
{
  Comp* aRootComp(theComps[0]);
  switch(LatticeType)
    {
    case HCP_LATTICE: 
      theAdjoiningVoxelSize = 12;
      theHCPk = theNormalizedVoxelRadius/sqrt(3); 
      theHCPh = theNormalizedVoxelRadius*sqrt(8.0/3);
      theHCPl = theNormalizedVoxelRadius*sqrt(3);
      break;
    case CUBIC_LATTICE:
      theAdjoiningVoxelSize = 6;
      break;
    }
  if(aRootComp->geometry == CUBOID)
    {
      //We do not give any leeway space between the simulation boundary
      //and the cell boundary if it is CUBOID to support
      //periodic boundary conditions:
      theCenterPoint.z = aRootComp->lengthZ/2; //row
      theCenterPoint.y = aRootComp->lengthY/2; //layer
      theCenterPoint.x = aRootComp->lengthX/2; //column
    }
  else
    {
      switch(LatticeType)
        {
        case HCP_LATTICE: 
          theCenterPoint.z = aRootComp->lengthZ/2+4*
            theNormalizedVoxelRadius; //row
          theCenterPoint.y = aRootComp->lengthY/2+2*theHCPl; //layer
          theCenterPoint.x = aRootComp->lengthX/2+2*theHCPh; //column
          break;
        case CUBIC_LATTICE:
          theCenterPoint.z = aRootComp->lengthZ/2+8*theNormalizedVoxelRadius;
          theCenterPoint.y = aRootComp->lengthY/2+8*theNormalizedVoxelRadius;
          theCenterPoint.x = aRootComp->lengthX/2+8*theNormalizedVoxelRadius;
          break;
        }
    }
  aRootComp->centerPoint = theCenterPoint; 
  switch(LatticeType)
    {
    case HCP_LATTICE: 
      theRowSize = (unsigned int)rint((theCenterPoint.z)/
                                      (theNormalizedVoxelRadius));
      theLayerSize = (unsigned int)rint((theCenterPoint.y*2)/theHCPl);
      theColSize = (unsigned int)rint((theCenterPoint.x*2)/theHCPh);
      break;
    case CUBIC_LATTICE:
      theRowSize = (unsigned int)rint((theCenterPoint.z)/
                                      (theNormalizedVoxelRadius));
      theLayerSize = (unsigned int)rint((theCenterPoint.y)/
                                      (theNormalizedVoxelRadius));
      theColSize = (unsigned int)rint((theCenterPoint.x)/
                                      (theNormalizedVoxelRadius));
      break;
    }
  //For the CUBOID cell geometrys, we need to readjust the size of
  //row, layer and column according to the boundary condition of its surfaces
  //to reflect the correct volume. This is because periodic boundary will
  //[consume a layer of the surface voxels:
  if(aRootComp->geometry == CUBOID)
    {
      //We need to increase the row, layer and col size by 2 because
      //the entire volume must be surrounded by nullID voxels to avoid
      //self-homodimerization reaction.
      theRowSize += 2;
      theColSize += 2;
      theLayerSize += 2;
      readjustSurfaceBoundarySizes(); 
    }

  //We use a single positive integer number for each voxel as its
  //coordinate.
  //Since we might use large number of voxels and the cell may grow in all
  //directions, it is best to use a starting coordinate that is half of
  //the largest positive integer supported:
  theStartCoord = UINT_MAX/2;
  //Drop back the start coordinate to the first voxel of the lattice:
  theStartCoord -= theStartCoord%(theRowSize*theLayerSize*
                                (theStartCoord/(theRowSize*theLayerSize)));
  theLattice.resize(theRowSize*theLayerSize*theColSize);
}

void SpatiocyteStepper::storeSimulationParameters()
{
  for(unsigned int i(0); i != theComps.size(); ++i)
    {
      Comp* aComp(theComps[i]); 
      if(aComp->dimension == 2)
        {
          aComp->actualArea =  (72*pow(VoxelRadius,2))*
            aComp->coords.size()/(6*pow(2,0.5)+4*pow(3,0.5)+
                                         3*pow(6, 0.5));
        }
      else // (aComp->dimension == 3)
        { 
          int voxelCnt(aComp->coords.size());
          for(unsigned int j(0); j != aComp->allSubs.size(); ++j)
            {
              voxelCnt += aComp->allSubs[j]->coords.size();
            }
          aComp->actualVolume = (4*pow(2,0.5)*pow(VoxelRadius,3))*
            voxelCnt;
        }
    }
}

void SpatiocyteStepper::printSimulationParameters()
{
  std::cout << std::endl;
  std::cout << "   Voxel radius, r_v:" << VoxelRadius << " m" << std::endl;
  std::cout << "   Simulation height:" << theCenterPoint.y*2*VoxelRadius*2 <<
    " m" << std::endl;
  std::cout << "   Simulation width:" << theCenterPoint.z*2*VoxelRadius*2 << 
    " m" << std::endl;
  std::cout << "   Simulation length:" << theCenterPoint.x*2*VoxelRadius*2 <<
    " m" << std::endl;
  std::cout << "   Layer size:" << theLayerSize << std::endl;
  std::cout << "   Row size:" << theRowSize << std::endl;
  std::cout << "   Column size:" << theColSize << std::endl;
  std::cout << "   Total allocated voxels:" << 
    theRowSize*theLayerSize*theColSize << std::endl;
  for(unsigned int i(0); i != theComps.size(); ++i)
    {
      Comp* aComp(theComps[i]);
      double aSpecVolume(aComp->specVolume);
      double aSpecArea(aComp->specArea);
      double anActualVolume(aComp->actualVolume);
      double anActualArea(aComp->actualArea);
      switch(aComp->geometry)
        {
        case CUBOID:
          std::cout << "   Cuboid ";
          break;
        case ELLIPSOID:
          std::cout << "   Ellipsoid ";
          break;
        case CYLINDER:
          std::cout << "   Cylinder (radius=" << aComp->lengthY*VoxelRadius
            << "m, length=" << (aComp->lengthX)*VoxelRadius*2 << "m) ";
          break;
        case ROD:
          std::cout << "   Rod (radius=" << aComp->lengthY*VoxelRadius << 
            "m, cylinder length=" <<
            (aComp->lengthX-aComp->lengthY*2)*VoxelRadius*2 << "m) ";
          break;
        }
      std::cout << aComp->system->getFullID().asString();
      switch (aComp->dimension)
      { 
      case 1:
          std::cout << " Line compartment:" << std::endl;
          break;
      case 2:
          std::cout << " Surface compartment:" << std::endl;
          std::cout << "     [" << int(aSpecArea*(6*sqrt(2)+4*sqrt(3)+3*sqrt(6))/
                              (72*VoxelRadius*VoxelRadius)) << 
            "] Specified surface voxels {n_s = S_specified*"
            << "(6*2^0.5+4*3^0.5+3*6^0.5)/(72*r_v^2}" << std::endl;
          std::cout << "     [" << aComp->coords.size() <<
            "] Actual surface voxels {n_s}" << std::endl;
          std::cout << "     [" << aSpecArea << " m^2] Specified surface area " <<
            "{S_specified}" << std::endl;
          std::cout << "     [" << anActualArea << " m^2] Actual surface area " <<
            "{S = (72*r_v^2)*n_s/(6*2^0.5+4*3^0.5+3*6^0.5)}" << std::endl;
          break;
      case 3:
      default:
          std::cout << " Volume compartment:" << std::endl;
          int voxelCnt(aComp->coords.size());
          for(unsigned int j(0); j != aComp->allSubs.size(); ++j)
            {
              voxelCnt += aComp->allSubs[j]->coords.size();
            }
          std::cout << "     [" << int(aSpecVolume/(4*sqrt(2)*pow(VoxelRadius, 3))) << 
            "] Specified volume voxels {n_v = V_specified/(4*2^0.5*r_v^3)}" <<
          std::endl;  
          std::cout << "     [" << voxelCnt << "] Actual volume voxels {n_v}"  << std::endl;
          std::cout << "     [" << aSpecVolume << " m^3] Specified volume {V_specified}"
            << std::endl; 
          std::cout << "     [" << anActualVolume << " m^3] Actual volume " <<
            "{V = (4*2^0.5*r_v^3)*n_v}" << std::endl; 
      }
    }
  std::cout << std::endl;
}

void SpatiocyteStepper::readjustSurfaceBoundarySizes()
{
  Comp* aRootComp(theComps[0]);
  //[XY, XZ, YZ]PLANE: the boundary type of the surface when 
  //the geometry of the root Comp is CUBOID.
  //Where the root cuboid compartment is enclosed with a surface compartment,
  //the boundary type can be either REFLECTIVE, REMOVE_UPPER, REMOVE_LOWER or
  //REMOVE_BOTH. To make the actualVolume of the root compartment equivalent
  //to the specVolume we need to increase the size of the row, col and layer
  //according to the additional voxels required to occupy the surface voxels.
  if(aRootComp->surfaceSub)
    {
      if(aRootComp->xyPlane == REFLECTIVE)
        {
          theRowSize += 2;
        }
      else if(aRootComp->xyPlane == REMOVE_UPPER ||
              aRootComp->xyPlane == REMOVE_LOWER)
        {
          theRowSize += 1;
        }
      if(aRootComp->yzPlane == REFLECTIVE)
        {
          theColSize += 2;
        }
      else if(aRootComp->yzPlane == REMOVE_UPPER ||
              aRootComp->yzPlane == REMOVE_LOWER)
        {
          theColSize += 1;
        }
      if(aRootComp->xzPlane == REFLECTIVE)
        {
          theLayerSize += 2;
        }
      else if(aRootComp->xzPlane == REMOVE_UPPER ||
              aRootComp->xzPlane == REMOVE_LOWER)
        {
          theLayerSize += 1;
        }
    }
  else
    {
  //Boundary type can also be either PERIODIC or REFLECTIVE when there is
  //no surface compartment for the root compartment.
  //Increase the size of [row,layer,col] by one voxel and make them odd sized
  //if the system uses periodic boundary conditions.
      if(aRootComp->yzPlane == PERIODIC)
        { 
          if(theColSize%2 != 1)
            {
              theColSize += 1;
            }
          else
            {
              theColSize += 2;
            }
        }
      if(aRootComp->xzPlane == PERIODIC)
        {
          if(theLayerSize%2 != 1)
            {
              theLayerSize +=1;
            }
          else
            {
              theLayerSize += 2;
            }
        }
      if(aRootComp->xyPlane == PERIODIC)
        {
          if(theRowSize%2 != 1)
            {
              theRowSize += 1;
            }
          else
            {
              theRowSize += 2;
            }
        }
      if(isPeriodicEdge)
        {
          if(theColSize%2 == 1)
            {
              theColSize += 1;
            }
          if(theLayerSize%2 == 1)
            {
              theLayerSize +=1;
            }
          if(theRowSize%2 == 1)
            {
              theRowSize += 1;
            }
        }
    }
}

void SpatiocyteStepper::constructLattice()
{ 
  Comp* aRootComp(theComps[0]);
  unsigned int aSize(theRowSize*theLayerSize*theColSize);
  unsigned int a(0);
  unsigned int b(theStartCoord);
  unsigned short rootID(aRootComp->vacantID);
  for(std::vector<Voxel>::iterator i(theLattice.begin()); a != aSize; ++i, ++a, ++b)
    { 
      (*i).adjoiningVoxels = new Voxel*[theAdjoiningVoxelSize];
      unsigned int aCol(a/(theRowSize*theLayerSize)); 
      unsigned int aLayer((a%(theRowSize*theLayerSize))/theRowSize); 
      unsigned int aRow((a%(theRowSize*theLayerSize))%theRowSize); 
      (*i).coord = b; 
      if(aRootComp->geometry == CUBOID || isInsideCoord(b, aRootComp, 0))
        {
          //By default, the voxel is vacant and we set it to the root id:
          (*i).id = rootID;
          for(unsigned int j(0); j != theAdjoiningVoxelSize; ++j)
            { 
              // By default let the adjoining voxel pointer point to the 
              // source voxel (i.e., itself)
              (*i).adjoiningVoxels[j] = &(*i);
            } 
          concatenateVoxel(&(*i), aRow, aLayer, aCol);
        }
      else
        {
          //We set id = theNullID if it is an invalid voxel, i.e., no molecules
          //will occupy it:
          (*i).id = theNullID;
          //Concatenate some of the null voxels close to the surface:
          if(isInsideCoord(b, aRootComp, 4))
            {
              concatenateVoxel(&(*i), aRow, aLayer, aCol);
            }
        }
    }
  if(aRootComp->geometry == CUBOID)
    {
      concatenatePeriodicSurfaces();
    }
}

void SpatiocyteStepper::setPeriodicEdge()
{
  isPeriodicEdge = true;
}

bool SpatiocyteStepper::isPeriodicEdgeCoord(unsigned int aCoord, Comp* aComp)
{
  unsigned int aRow;
  unsigned int aLayer;
  unsigned int aCol;
  coord2global(aCoord, &aRow, &aLayer, &aCol);
  if(aComp->system->getSuperSystem()->isRootSystem() &&
     ((aRow <= 1 && aCol <= 1) ||
      (aRow <= 1 && aLayer <= 1) ||
      (aCol <= 1 && aLayer <= 1) ||
      (aRow <= 1 && aCol >= theColSize-2) ||
      (aRow <= 1 && aLayer >= theLayerSize-2) ||
      (aCol <= 1 && aRow >= theRowSize-2) ||
      (aCol <= 1 && aLayer >= theLayerSize-2) ||
      (aLayer <= 1 && aRow >= theRowSize-2) ||
      (aLayer <= 1 && aCol >= theColSize-2) ||
      (aRow >= theRowSize-2 && aCol >= theColSize-2) ||
      (aRow >= theRowSize-2 && aLayer >= theLayerSize-2) ||
      (aCol >= theColSize-2 && aLayer >= theLayerSize-2)))
    {
      return true;
    }
  return false;
}

bool SpatiocyteStepper::isRemovableEdgeCoord(unsigned int aCoord, Comp* aComp)
{
  unsigned int aRow;
  unsigned int aLayer;
  unsigned int aCol;
  coord2global(aCoord, &aRow, &aLayer, &aCol);
  int sharedCnt(0);
  int removeCnt(0);
  //Minus 1 to maxRow to account for surfaces that use two rows to envelope the
  //volume:
  if(aRow >= aComp->maxRow-1)
    {
      ++sharedCnt;
      if(aComp->xyPlane == REMOVE_UPPER || aComp->xyPlane == REMOVE_BOTH)
        {
          ++removeCnt;
        }
    }
  //Add 1 to minRow to account for surfaces that use two rows to envelope the
  //volume:
  if(aRow <= aComp->minRow+1)
    {
      ++sharedCnt;
      if(aComp->xyPlane == REMOVE_LOWER || aComp->xyPlane == REMOVE_BOTH)
        {
          ++removeCnt;
        }
    }
  if(aLayer >= aComp->maxLayer-1)
    {
      ++sharedCnt;
      if(aComp->xzPlane == REMOVE_UPPER || aComp->xzPlane == REMOVE_BOTH)
        {
          ++removeCnt;
        }
    }
  if(aLayer <= aComp->minLayer+1)
    {
      ++sharedCnt;
      if(aComp->xzPlane == REMOVE_LOWER || aComp->xzPlane == REMOVE_BOTH)
        {
          ++removeCnt;
        }
    }
  if(aCol >= aComp->maxCol)
    {
      ++sharedCnt;
      if(aComp->yzPlane == REMOVE_UPPER || aComp->yzPlane == REMOVE_BOTH)
        {
          ++removeCnt;
        }
    }
  if(aCol <= aComp->minCol) 
    {
      ++sharedCnt;
      if(aComp->yzPlane == REMOVE_LOWER || aComp->yzPlane == REMOVE_BOTH)
        {
          ++removeCnt;
        }
    }
  if(!removeCnt)
    {
      return false;
    }
  else
    {
      return sharedCnt == removeCnt;
    }
}

void SpatiocyteStepper::concatenateVoxel(Voxel* aVoxel,
                                         unsigned int aRow,
                                         unsigned int aLayer,
                                         unsigned int aCol)
{
  if(aRow > 0)
    { 
      concatenateRows(aVoxel, aRow-1, aLayer, aCol); 
    } 
  if(aLayer > 0)
    {
      concatenateLayers(aVoxel, aRow, aLayer-1, aCol); 
    }
  if(aCol > 0)
    { 
      concatenateCols(aVoxel, aRow, aLayer, aCol-1); 
    }
}

Variable* SpatiocyteStepper::getVariable(System* aSystem, String const& anID)
{
  FOR_ALL(System::Variables, aSystem->getVariables())
    {
      Variable* aVariable(i->second);
      if(aVariable->getID() == anID)
        {
          return aVariable;
        }
    }
  return NULL;
}

void SpatiocyteStepper::setCompProperties(Comp* aComp)
{
  System* aSystem(aComp->system);
  double aRadius(0);
  double aLongRadius(0);
  switch(aComp->geometry)
    {
    case CUBOID:
      if(!aComp->lengthX)
        {
          THROW_EXCEPTION(NotFound,
                          "Property LENGTHX of the Cuboid Comp " +
                          aSystem->getFullID().asString() + " is not defined.");
        }
      if(!aComp->lengthY)
        {
          THROW_EXCEPTION(NotFound,
                          "Property LENGTHY of the Cuboid Comp " +
                          aSystem->getFullID().asString() + " is not defined.");
        }
      if(!aComp->lengthZ)
        {
          THROW_EXCEPTION(NotFound,
                          "Property LENGTHZ of the Cuboid Comp " +
                          aSystem->getFullID().asString() + " is not defined.");
        }
      aComp->specVolume = aComp->lengthX*aComp->lengthY*
        aComp->lengthZ;
      aComp->specArea = getCuboidSpecArea(aComp);
      break;
    case ELLIPSOID:
      if(!aComp->lengthX)
        {
          THROW_EXCEPTION(NotFound,
                          "Property LENGTHX of the Ellipsoid Comp " +
                          aSystem->getFullID().asString() + " is not defined.");
        }
      if(!aComp->lengthY)
        {
          THROW_EXCEPTION(NotFound,
                          "Property LENGTHY of the Ellipsoid Comp " +
                          aSystem->getFullID().asString() + " is not defined.");
        }
      if(!aComp->lengthZ)
        {
          THROW_EXCEPTION(NotFound,
                          "Property LENGTHZ of the Ellipsoid Comp " +
                          aSystem->getFullID().asString() + " is not defined.");
        }
      aComp->specVolume = 4*M_PI*aComp->lengthX*
        aComp->lengthY* aComp->lengthZ/24;
      aComp->specArea = 4*M_PI*
        pow((pow(aComp->lengthX/2, 1.6075)*
             pow(aComp->lengthY/2, 1.6075)+
             pow(aComp->lengthX/2, 1.6075)*
             pow(aComp->lengthZ/2, 1.6075)+
             pow(aComp->lengthY/2, 1.6075)*
             pow(aComp->lengthZ/2, 1.6075))/ 3, 1/1.6075); 
      break;
    case CYLINDER:
      if(!aComp->lengthX)
        {
          THROW_EXCEPTION(NotFound, "Property LENGTHX of the Cylinder Comp "
                          + aSystem->getFullID().asString() + " not defined." );
        }
      if(!aComp->lengthY)
        {
          THROW_EXCEPTION(NotFound, "Property LENGTHY of the Cylinder Comp "
                          + aSystem->getFullID().asString() + ", which "
                          + "specifies its diameter, not defined." );
        }
      aRadius = aComp->lengthY/2;
      aComp->specVolume = (aComp->lengthX)*(M_PI*aRadius*aRadius);
      aComp->lengthZ = aComp->lengthY;
      aComp->specArea = 2*M_PI*aRadius*(aComp->lengthX);
      break;
    case ROD:
      if(!aComp->lengthX)
        {
          THROW_EXCEPTION(NotFound, "Property LENGTHX of the Rod Comp "
                          + aSystem->getFullID().asString() + ", which "
                          + "specifies the rod length (cylinder length + "
                          + "2*hemisphere radius), not defined." );
        }
      if(!aComp->lengthY)
        {
          THROW_EXCEPTION(NotFound, "Property LENGTHY of the Rod Comp "
                          + aSystem->getFullID().asString() + ", which "
                          + "specifies its diameter, not defined." );
        }
      aRadius = aComp->lengthY/2;
      aComp->specVolume = (aComp->lengthX+(4*aRadius/3)-(2*aRadius))*
        (M_PI*aRadius*aRadius);
      aComp->lengthZ = aComp->lengthY;
      aComp->specArea = 4*M_PI*aRadius*aRadius+
        2*M_PI*aRadius*(aComp->lengthX-2*aRadius);
      break;
    case TORUS:
      if(!aComp->lengthX)
        {
          THROW_EXCEPTION(NotFound, "Property LENGTHX of the Torus Comp "
                          + aSystem->getFullID().asString() + ", which "
                          + "specifies the diameter of the torus, "
                          + "not defined." );
        }
      if(!aComp->lengthY)
        {
          THROW_EXCEPTION(NotFound, "Property LENGTHY of the Torus Comp "
                          + aSystem->getFullID().asString() + ", which "
                          + "specifies the diameter of tube, not defined." );
        }
      aRadius = aComp->lengthY/2;
      aLongRadius = (aComp->lengthX-aComp->lengthY)/2;
      aComp->lengthZ = aComp->lengthX;
      aComp->specVolume = M_PI*aRadius*aRadius*2*M_PI*aLongRadius;
      aComp->specArea = 2*M_PI*aRadius*2*M_PI*aLongRadius;
      break;
    case PYRAMID:
      if(!aComp->lengthX)
        {
          THROW_EXCEPTION(NotFound, "Property LENGTHX of the Pyramid Comp "
                          + aSystem->getFullID().asString() + ", which "
                          + "specifies the length of the pyramid, "
                          + "not defined." );
        }
      if(!aComp->lengthY)
        {
          THROW_EXCEPTION(NotFound, "Property LENGTHY of the Pyramid Comp "
                          + aSystem->getFullID().asString() + ", which "
                          + "specifies the height of the pyramid, "
                          + "not defined." );
        }
      if(!aComp->lengthZ)
        {
          THROW_EXCEPTION(NotFound, "Property LENGTHZ of the Pyramid Comp "
                          + aSystem->getFullID().asString() + ", which "
                          + "specifies the depth of the pyramid, "
                          + "not defined." );
        }
      aComp->specVolume = aComp->lengthX*aComp->lengthY*aComp->lengthZ/3;
      aRadius = sqrt(pow(aComp->lengthX/2,2)+pow(aComp->lengthY,2));
      aLongRadius = sqrt(pow(aComp->lengthZ/2,2)+pow(aComp->lengthY,2));
      aComp->specArea =  aComp->lengthZ*aComp->lengthX+
        2*(aComp->lengthZ*aRadius)+2*(aComp->lengthX*aRadius);
      break;
    }
  aComp->lengthX /= VoxelRadius*2;
  aComp->lengthY /= VoxelRadius*2;
  aComp->lengthZ /= VoxelRadius*2;
}

void SpatiocyteStepper::setCompCenterPoint(Comp* aComp)
{
  System* aSystem(aComp->system);
  System* aSuperSystem(aSystem->getSuperSystem());
  if(aComp->dimension == 2)
    {
      aSystem = aComp->system->getSuperSystem();
      aSuperSystem = aSystem->getSuperSystem();
    }
  else if(aComp->dimension == 1)
    {
      aSystem = aComp->system->getSuperSystem()->getSuperSystem();
      aSuperSystem = aSystem->getSuperSystem();
    }
  Comp* aSuperComp(system2Comp(aSuperSystem));
  //The center with reference to the immediate super system:
  aComp->centerPoint = aSuperComp->centerPoint;
  aComp->centerPoint.x += aComp->originX*aSuperComp->lengthX/2;
  aComp->centerPoint.y += aComp->originY*aSuperComp->lengthY/2;
  aComp->centerPoint.z += aComp->originZ*aSuperComp->lengthZ/2;
}

double SpatiocyteStepper::getCuboidSpecArea(Comp* aComp)
{
  double anArea(0);
  if(aComp->xyPlane == UNIPERIODIC || 
     aComp->xyPlane == REMOVE_UPPER ||
     aComp->xyPlane == REMOVE_LOWER)
    { 
      anArea += aComp->lengthX*aComp->lengthY; 
    }
  else if(aComp->xyPlane == REFLECTIVE)
    {
      anArea += 2*aComp->lengthX*aComp->lengthY; 
    }
  if(aComp->xzPlane == UNIPERIODIC || 
     aComp->xzPlane == REMOVE_UPPER ||
     aComp->xzPlane == REMOVE_LOWER)
    { 
      anArea += aComp->lengthX*aComp->lengthZ; 
    }
  else if(aComp->xzPlane == REFLECTIVE)
    {
      anArea += 2*aComp->lengthX*aComp->lengthZ; 
    }
  if(aComp->yzPlane == UNIPERIODIC || 
     aComp->yzPlane == REMOVE_UPPER ||
     aComp->yzPlane == REMOVE_LOWER)
    { 
      anArea += aComp->lengthY*aComp->lengthZ; 
    }
  else if(aComp->yzPlane == REFLECTIVE)
    {
      anArea += 2*aComp->lengthY*aComp->lengthZ; 
    }
  return anArea;;
}

Point SpatiocyteStepper::coord2point(unsigned int aCoord)
{
  unsigned int aGlobalCol;
  unsigned int aGlobalLayer;
  unsigned int aGlobalRow;
  coord2global(aCoord, &aGlobalRow, &aGlobalLayer, &aGlobalCol);
  //the center point of a voxel 
  Point aPoint;
  switch(LatticeType)
    {
    case HCP_LATTICE: 
      aPoint.y = (aGlobalCol%2)*theHCPk+theHCPl*aGlobalLayer;
      aPoint.z = aGlobalRow*2*theNormalizedVoxelRadius+
        ((aGlobalLayer+aGlobalCol)%2)*theNormalizedVoxelRadius;
      aPoint.x = aGlobalCol*theHCPh;
      break;
    case CUBIC_LATTICE:
      aPoint.y = aGlobalLayer*2*theNormalizedVoxelRadius;
      aPoint.z = aGlobalRow*2*theNormalizedVoxelRadius;
      aPoint.x = aGlobalCol*2*theNormalizedVoxelRadius;
      break;
    }
  return aPoint;
}

Voxel* SpatiocyteStepper::point2voxel(Point aPoint)
{
  unsigned int aGlobalCol(0);
  unsigned int aGlobalLayer(0);
  unsigned int aGlobalRow(0);
  switch(LatticeType)
    {
    case HCP_LATTICE: 
      aGlobalCol = (unsigned int)(aPoint.x/theHCPh);
      aGlobalLayer = (unsigned int)((aPoint.y-(aGlobalCol%2)*theHCPk)/theHCPl);
      aGlobalRow = (unsigned int)((aPoint.z-((aGlobalLayer+aGlobalCol)%2)*
          theNormalizedVoxelRadius)/(2*theNormalizedVoxelRadius));
      break;
    case CUBIC_LATTICE:
      aGlobalCol = (unsigned int)(aPoint.x/(2*theNormalizedVoxelRadius));
      aGlobalLayer = (unsigned int)(aPoint.y/(2*theNormalizedVoxelRadius));
      aGlobalRow = (unsigned int)(aPoint.z/(2*theNormalizedVoxelRadius));
      break;
    }
  if(aGlobalCol < 0)
    {
      aGlobalCol = 0;
    }
  else if(aGlobalCol >= theColSize)
    {
      aGlobalCol = theColSize-1;
    }
  if(aGlobalRow < 0)
    {
      aGlobalRow = 0;
    }
  else if(aGlobalRow >= theRowSize)
    {
      aGlobalRow = theRowSize-1;
    }
  if(aGlobalLayer < 0)
    {
      aGlobalLayer = 0;
    }
  else if(aGlobalLayer >= theLayerSize)
    {
      aGlobalLayer = theLayerSize-1;
    }
  return &theLattice[aGlobalRow+
                     theRowSize*aGlobalLayer+
                     theRowSize*theLayerSize*aGlobalCol];
}

void SpatiocyteStepper::coord2global(unsigned int aCoord,
                                     unsigned int* aGlobalRow,
                                     unsigned int* aGlobalLayer,
                                     unsigned int* aGlobalCol) 
{
  *aGlobalCol = (aCoord-theStartCoord)/(theRowSize*theLayerSize);
  *aGlobalLayer = ((aCoord-theStartCoord)%(theRowSize*theLayerSize))/theRowSize;
  *aGlobalRow = ((aCoord-theStartCoord)%(theRowSize*theLayerSize))%theRowSize;
}

void SpatiocyteStepper::concatenateLayers(Voxel* aVoxel,
                                          unsigned int aRow,
                                          unsigned int aLayer,
                                          unsigned int aCol)
{
  Voxel* anAdjoiningVoxel(&theLattice[
                          aRow+
                          theRowSize*aLayer+
                          theRowSize*theLayerSize*aCol]);
  switch(LatticeType)
    {
    case HCP_LATTICE: 
      if((aLayer+1)%2+(aCol)%2 == 1)
        {
          aVoxel->adjoiningVoxels[VENTRALN] = anAdjoiningVoxel;
          anAdjoiningVoxel->adjoiningVoxels[DORSALS] = aVoxel;
          if(aRow < theRowSize-1)
            {
              anAdjoiningVoxel = (&theLattice[
                                  aRow+1+ 
                                  theRowSize*aLayer+
                                  theRowSize*theLayerSize*aCol]);
              aVoxel->adjoiningVoxels[VENTRALS] = anAdjoiningVoxel;
              anAdjoiningVoxel->adjoiningVoxels[DORSALN] = aVoxel;
            }
        }
      else
        {
          aVoxel->adjoiningVoxels[VENTRALS] = anAdjoiningVoxel;
          anAdjoiningVoxel->adjoiningVoxels[DORSALN] = aVoxel;
          if(aRow > 0)
            {
              anAdjoiningVoxel = (&theLattice[
                                  aRow-1+ 
                                  theRowSize*aLayer+
                                  theRowSize*theLayerSize*aCol]);
              aVoxel->adjoiningVoxels[VENTRALN] = anAdjoiningVoxel;
              anAdjoiningVoxel->adjoiningVoxels[DORSALS] = aVoxel;
            }
        }
      break;
    case CUBIC_LATTICE: 
      aVoxel->adjoiningVoxels[VENTRAL] = anAdjoiningVoxel;
      anAdjoiningVoxel->adjoiningVoxels[DORSAL] = aVoxel;
      break;
    }
}


void SpatiocyteStepper::concatenateCols(Voxel* aVoxel,
                                        unsigned int aRow,
                                        unsigned int aLayer,
                                        unsigned int aCol)
{
  Voxel* anAdjoiningVoxel(&theLattice[
                          aRow+
                          theRowSize*aLayer+
                          theRowSize*theLayerSize*aCol]);
  switch(LatticeType)
    {
    case HCP_LATTICE: 
      if(aLayer%2 == 0)
        {
          if((aCol+1)%2 == 1)
            {
              aVoxel->adjoiningVoxels[NW] = anAdjoiningVoxel;
              anAdjoiningVoxel->adjoiningVoxels[SE] = aVoxel;
              if(aRow < theRowSize - 1)
                {
                  anAdjoiningVoxel = (&theLattice[
                                      aRow+1+ 
                                      theRowSize*aLayer+
                                      theRowSize*theLayerSize*aCol]);
                  aVoxel->adjoiningVoxels[SW] = anAdjoiningVoxel;
                  anAdjoiningVoxel->adjoiningVoxels[NE] = aVoxel;
                }
              if(aLayer < theLayerSize-1)
                {
                  anAdjoiningVoxel = (&theLattice[
                                      aRow+ 
                                      theRowSize*(aLayer+1)+
                                      theRowSize*theLayerSize*aCol]);
                  aVoxel->adjoiningVoxels[WEST] = anAdjoiningVoxel;
                  anAdjoiningVoxel->adjoiningVoxels[EAST] = aVoxel;
                }
            }
          else
            {
              aVoxel->adjoiningVoxels[SW] = anAdjoiningVoxel;
              anAdjoiningVoxel->adjoiningVoxels[NE] = aVoxel;
              if(aRow > 0)
                {
                  anAdjoiningVoxel = (&theLattice[
                                      aRow-1+ 
                                      theRowSize*aLayer+
                                      theRowSize*theLayerSize*aCol]);
                  aVoxel->adjoiningVoxels[NW] = anAdjoiningVoxel;
                  anAdjoiningVoxel->adjoiningVoxels[SE] = aVoxel;
                }
              if(aLayer > 0)
                {
                  anAdjoiningVoxel = (&theLattice[
                                      aRow+ 
                                      theRowSize*(aLayer-1)+
                                      theRowSize*theLayerSize*aCol]);
                  aVoxel->adjoiningVoxels[WEST] = anAdjoiningVoxel;
                  anAdjoiningVoxel->adjoiningVoxels[EAST] = aVoxel;
                }
            }
        }
      else
        {
          if((aCol+1)%2 == 1)
            {
              aVoxel->adjoiningVoxels[SW] = anAdjoiningVoxel;
              anAdjoiningVoxel->adjoiningVoxels[NE] = aVoxel;
              if(aRow > 0)
                {
                  anAdjoiningVoxel = (&theLattice[
                                      aRow-1+ 
                                      theRowSize*aLayer+
                                      theRowSize*theLayerSize*aCol]);
                  aVoxel->adjoiningVoxels[NW] = anAdjoiningVoxel;
                  anAdjoiningVoxel->adjoiningVoxels[SE] = aVoxel;
                }
              if(aLayer < theLayerSize-1)
                {
                  anAdjoiningVoxel = (&theLattice[
                                      aRow+ 
                                      theRowSize*(aLayer+1)+
                                      theRowSize*theLayerSize*aCol]);
                  aVoxel->adjoiningVoxels[WEST] = anAdjoiningVoxel;
                  anAdjoiningVoxel->adjoiningVoxels[EAST] = aVoxel;
                }
            }
          else
            {
              aVoxel->adjoiningVoxels[NW] = anAdjoiningVoxel;
              anAdjoiningVoxel->adjoiningVoxels[SE] = aVoxel;
              if(aRow < theRowSize - 1)
                {
                  anAdjoiningVoxel = (&theLattice[
                                      aRow+1+ 
                                      theRowSize*aLayer+
                                      theRowSize*theLayerSize*aCol]);
                  aVoxel->adjoiningVoxels[SW] = anAdjoiningVoxel;
                  anAdjoiningVoxel->adjoiningVoxels[NE] = aVoxel;
                }
              if(aLayer > 0)
                {
                  anAdjoiningVoxel = (&theLattice[
                                      aRow+ 
                                      theRowSize*(aLayer-1)+
                                      theRowSize*theLayerSize*aCol]);
                  aVoxel->adjoiningVoxels[WEST] = anAdjoiningVoxel;
                  anAdjoiningVoxel->adjoiningVoxels[EAST] = aVoxel;
                }
            }
        }
      break;
    case CUBIC_LATTICE: 
      aVoxel->adjoiningVoxels[WEST] = anAdjoiningVoxel;
      anAdjoiningVoxel->adjoiningVoxels[EAST] = aVoxel;
      break;
    }
}

void SpatiocyteStepper::concatenateRows(Voxel* aVoxel,
                                        unsigned int aRow,
                                        unsigned int aLayer,
                                        unsigned int aCol)
{
  Voxel* anAdjoiningVoxel(&theLattice[
                          aRow+ 
                          theRowSize*aLayer+ 
                          theRowSize*theLayerSize*aCol]);
  aVoxel->adjoiningVoxels[NORTH] = anAdjoiningVoxel;
  anAdjoiningVoxel->adjoiningVoxels[SOUTH] = aVoxel;
}


void SpatiocyteStepper::concatenatePeriodicSurfaces()
{
  Comp* aRootComp(theComps[0]);
  for(unsigned int i(0); i<=theRowSize*theLayerSize*theColSize-theRowSize;
      i+=theRowSize)
    {
      unsigned int j(i+theRowSize-1);
      Voxel* aSrcVoxel(&theLattice[
                       coord2row(i)+ 
                       theRowSize*coord2layer(i)+ 
                       theRowSize*theLayerSize*coord2col(i)]); 
      Voxel* aDestVoxel(&theLattice[
                        coord2row(j)+
                        theRowSize*coord2layer(i)+
                        theRowSize*theLayerSize*coord2col(i)]); 
      if(aRootComp->xyPlane == UNIPERIODIC)
        { 
          replaceUniVoxel(aSrcVoxel, aDestVoxel);
        }
      else if(aRootComp->xyPlane == PERIODIC)
        { 
          replaceVoxel(aSrcVoxel, aDestVoxel);
        }
      else if(!isPeriodicEdge)
        {
          //We cannot have valid voxels pointing to itself it it is not periodic
          //to avoid incorrect homodimerization reaction. So we set such
          //molecules to null ID.
          aSrcVoxel->id = theNullID;
          aDestVoxel->id = theNullID;
        }
    }
  for(unsigned int i(0); i<=theRowSize*theLayerSize*(theColSize-1)+theRowSize;)
    {
      unsigned int j(theRowSize*(theLayerSize-1)+i);
      Voxel* aSrcVoxel(&theLattice[
                       coord2row(i)+
                       theRowSize*coord2layer(i)+ 
                       theRowSize*theLayerSize*coord2col(i)]); 
      Voxel* aDestVoxel(&theLattice[
                        coord2row(i)+
                        theRowSize*coord2layer(j)+ 
                        theRowSize*theLayerSize*coord2col(i)]); 

      if(aRootComp->xzPlane == UNIPERIODIC)
        { 
          replaceUniVoxel(aSrcVoxel, aDestVoxel);
        }
      else if(aRootComp->xzPlane == PERIODIC)
        { 
          replaceVoxel(aSrcVoxel, aDestVoxel);
        }
      else if(!isPeriodicEdge)
        {
          aSrcVoxel->id = theNullID;
          aDestVoxel->id = theNullID;
        }
      ++i;
      if(coord2layer(i) != 0)
        {
          i += (theLayerSize-1)*theRowSize;
        }
    }
  for(unsigned int i(0); i!=theRowSize*theLayerSize; ++i)
    {
      Voxel* aSrcVoxel(&theLattice[
                       coord2row(i)+
                       theRowSize*coord2layer(i)+ 
                       theRowSize*theLayerSize*0]); 
      Voxel* aDestVoxel(&theLattice[
                        coord2row(i)+
                        theRowSize*coord2layer(i)+ 
                        theRowSize*theLayerSize*(theColSize-1)]); 
      if(aRootComp->yzPlane == UNIPERIODIC)
        { 
          replaceUniVoxel(aSrcVoxel, aDestVoxel);
        }
      else if(aRootComp->yzPlane == PERIODIC)
        { 
          replaceVoxel(aSrcVoxel, aDestVoxel);
        }
      else if(!isPeriodicEdge)
        {
          aSrcVoxel->id = theNullID;
          aDestVoxel->id = theNullID;
        }
    }
}

unsigned int SpatiocyteStepper::coord2row(unsigned int aCoord)
{
  return (aCoord%(theRowSize*theLayerSize))%theRowSize;
}

unsigned int SpatiocyteStepper::coord2layer(unsigned int aCoord)
{
  return (aCoord%(theRowSize*theLayerSize))/theRowSize;
}

unsigned int SpatiocyteStepper::coord2col(unsigned int aCoord)
{
  return aCoord/(theRowSize*theLayerSize);
}

void SpatiocyteStepper::replaceVoxel(Voxel* aSrcVoxel, Voxel* aDestVoxel)
{
  if(aSrcVoxel->id != theNullID && aDestVoxel->id != theNullID)
    {
      for(unsigned int j(0); j!=theAdjoiningVoxelSize; ++j)
        {
          if(aSrcVoxel->adjoiningVoxels[j] == aSrcVoxel &&
              aDestVoxel->adjoiningVoxels[j] != aDestVoxel)
            {
              aSrcVoxel->adjoiningVoxels[j] = aDestVoxel->adjoiningVoxels[j];
              for(unsigned int k(0); k!=theAdjoiningVoxelSize; ++k)
                {
                  if(aDestVoxel->adjoiningVoxels[j]->adjoiningVoxels[k] ==
                     aDestVoxel)
                    {
                      aDestVoxel->adjoiningVoxels[j]->adjoiningVoxels[k] =
                        aSrcVoxel;
                    }
                }
            }
        }
      aDestVoxel->id = theNullID;
    }
}

void SpatiocyteStepper::replaceUniVoxel(Voxel* aSrcVoxel, Voxel* aDestVoxel)
{
  if(aSrcVoxel->id != theNullID && aDestVoxel->id != theNullID)
    {
      for(unsigned int j(0); j!=theAdjoiningVoxelSize; ++j)
        {
          if(aSrcVoxel->adjoiningVoxels[j] == aSrcVoxel)
            {
              for(unsigned int k(0); k!=theAdjoiningVoxelSize; ++k)
                {
                  if(aDestVoxel->adjoiningVoxels[j]->adjoiningVoxels[k] ==
                     aDestVoxel)
                    {
                      aDestVoxel->adjoiningVoxels[j]->adjoiningVoxels[k] =
                        aSrcVoxel;
                    }
                }
            }
        }
      aDestVoxel->id = theNullID;
    }
}

void SpatiocyteStepper::shuffleAdjoiningVoxels()
{
  for(std::vector<Voxel>::iterator i(theLattice.begin()); i != theLattice.end(); ++i)
    {
      if((*i).id != theNullID)
        { 
          gsl_ran_shuffle(getRng(), (*i).adjoiningVoxels, theAdjoiningVoxelSize,
                          sizeof(Voxel*));
        }
    }
}

void SpatiocyteStepper::setCompVoxelProperties()
{
  for(std::vector<Comp*>::iterator i(theComps.begin());
      i != theComps.end(); ++i)
    {
      switch ((*i)->dimension)
        {
        case 1:
          setLineCompProperties(*i);
          setLineVoxelProperties(*i);
          break;
        case 2:
          setSurfaceCompProperties(*i);
          setSurfaceVoxelProperties(*i);
          break;
        case 3:
        default:
            setVolumeCompProperties(*i);
        }
    }
}

void SpatiocyteStepper::setLineCompProperties(Comp* aComp)
{
    setSurfaceCompProperties(aComp);
}

void SpatiocyteStepper::setLineVoxelProperties(Comp* aComp)
{
    setSurfaceVoxelProperties(aComp);
}

void SpatiocyteStepper::setSurfaceCompProperties(Comp* aComp)
{
  removePeriodicEdgeVoxels(aComp);
  removeSurfaces(aComp);
  setDiffusiveComp(aComp);
}

void SpatiocyteStepper::setVolumeCompProperties(Comp* aComp)
{
  setDiffusiveComp(aComp);
}

void SpatiocyteStepper::setSurfaceVoxelProperties(Comp* aComp)
{
  if(!aComp->diffusiveComp)
    {
      for(std::vector<unsigned int>::iterator i(aComp->coords.begin());
          i != aComp->coords.end(); ++i)
        {
          Voxel* aVoxel(&theLattice[*i]);
          optimizeSurfaceVoxel(aVoxel, aComp);
          setSurfaceSubunit(aVoxel, aComp);
        }
    }
}

void SpatiocyteStepper::setDiffusiveComp(Comp* aComp)
{
  FOR_ALL(System::Variables, aComp->system->getVariables())
    {
      Variable* aVariable(i->second);
      if(aVariable->getID() == "DIFFUSIVE")
        {
          String aStringID(aVariable->getName()); 
          aStringID = "System:" + aStringID;
          FullID aFullID(aStringID);
          System* aSystem(static_cast<System*>(getModel()->getEntity(aFullID)));
          aComp->diffusiveComp = system2Comp(aSystem);
        }
    }
  if(aComp->diffusiveComp)
    {
      for(std::vector<unsigned int>::iterator j(aComp->coords.begin());
          j != aComp->coords.end(); ++j)
        {
          Voxel* aVoxel(&theLattice[*j]);
          aVoxel->id = aComp->diffusiveComp->vacantID;
          aComp->diffusiveComp->coords.push_back( aVoxel->coord-theStartCoord);
        }
      aComp->coords.clear();
    }
}


void SpatiocyteStepper::removePeriodicEdgeVoxels(Comp* aComp)
{ 
  if(isPeriodicEdge)
    {
      std::vector<unsigned int> coords;
      for(std::vector<unsigned int>::iterator j(aComp->coords.begin());
          j != aComp->coords.end(); ++j)
        {
          Voxel* aVoxel(&theLattice[*j]);
          if(isPeriodicEdgeCoord(aVoxel->coord, aComp))
            {
              aVoxel->id = theNullID;
            }
          else
            {
              coords.push_back(*j);
            }
        }
      aComp->coords = coords;
    }
}

void SpatiocyteStepper::removeSurfaces(Comp* aComp)
{ 
  std::vector<unsigned int> coords;
  for(std::vector<unsigned int>::iterator j(aComp->coords.begin());
      j != aComp->coords.end(); ++j)
    {
      Voxel* aVoxel(&theLattice[*j]);
      if(isRemovableEdgeCoord(aVoxel->coord, aComp))
        { 
          Comp* aSuperComp(system2Comp(aComp->system->getSuperSystem())); 
          aVoxel->id = aSuperComp->vacantID;
          aSuperComp->coords.push_back(*j);
        }
      else
        {
          coords.push_back(*j);
        }
    }
  aComp->coords = coords;
}

void SpatiocyteStepper::optimizeSurfaceVoxel(Voxel* aVoxel,
                                             Comp* aComp)
{
  unsigned short surfaceID(aComp->vacantID);
  aVoxel->surfaceVoxels = new std::vector<std::vector<Voxel*> >;
  aVoxel->surfaceVoxels->resize(4);
  std::vector<Voxel*>& immediateSurface((*aVoxel->surfaceVoxels)[IMMEDIATE]);
  std::vector<Voxel*>& extendedSurface((*aVoxel->surfaceVoxels)[EXTENDED]);
  std::vector<Voxel*>& innerVolume((*aVoxel->surfaceVoxels)[INNER]);
  std::vector<Voxel*>& outerVolume((*aVoxel->surfaceVoxels)[OUTER]);
  std::vector<std::vector<Voxel*> > sharedVoxelsList;
  Voxel** forward(aVoxel->adjoiningVoxels);
  Voxel** reverse(forward+theAdjoiningVoxelSize);
  std::vector<Voxel*> adjoiningCopy;
  for(unsigned int k(0); k != theAdjoiningVoxelSize; ++k)
    {
      adjoiningCopy.push_back(forward[k]);
    }
  //Separate adjoining surface voxels and adjoining volume voxels.
  //Put the adjoining surface voxels at the beginning of the
  //adjoiningVoxels list while the volume voxels are put at the end:
  for(std::vector<Voxel*>::iterator l(adjoiningCopy.begin());
      l != adjoiningCopy.end(); ++l)
    {
      if((*l) != aVoxel && (*l)->id != theNullID 
         && id2Comp((*l)->id)->dimension <= aComp->dimension)
        {
          (*forward) = (*l);
          ++forward;
          //immediateSurface contains all adjoining surface voxels except the 
          //source voxel, aVoxel:
          immediateSurface.push_back(*l);
          for(unsigned int m(0); m != theAdjoiningVoxelSize; ++m)
            {
              //extendedSurface contains the adjoining surface voxels of
              //adjoining surface voxels. They do not include the source voxel
              //and its adjoining voxels:
              Voxel* extendedVoxel((*l)->adjoiningVoxels[m]);
              if(extendedVoxel->id == surfaceID && extendedVoxel != aVoxel &&
                 std::find(adjoiningCopy.begin(), adjoiningCopy.end(),
                      extendedVoxel) == adjoiningCopy.end())
                {
                  std::vector<Voxel*>::iterator n(std::find(extendedSurface.begin(),
                        extendedSurface.end(), extendedVoxel));
                  if(n == extendedSurface.end())
                    {
                      extendedSurface.push_back(extendedVoxel);
                      //We require shared immediate voxel which
                      //connects the extended voxel with the source voxel 
                      //for polymerization. Create a new list of shared
                      //immediate voxel each time a new extended voxel is added:
                      std::vector<Voxel*> sharedVoxels;
                      sharedVoxels.push_back(*l);
                      sharedVoxelsList.push_back(sharedVoxels);
                    }
                  else
                    {
                      //An extended voxel may have multiple shared immediate
                      //voxels, so we insert the additional ones in the list:
                      sharedVoxelsList[n-extendedSurface.begin()].push_back(*l);
                    }
                }
            }
        }
      else
        {
          --reverse;
          (*reverse) = (*l);
          //We know that it is not a surface voxel, so it would not
          //be a self-pointed adjoining voxel. If it is not inside the
          //surface (i.e., the parent volume) Comp, it must be an
          //outer volume voxel:
          if(!isInsideCoord((*l)->coord, aComp, 0))
            {
              outerVolume.push_back(*l);
            }
          //otherwise, it must be an inner volume voxel:
          else
            {
              innerVolume.push_back(*l);
            }
        }
    } 
  for(std::vector<std::vector<Voxel*> >::iterator i(sharedVoxelsList.begin());
      i != sharedVoxelsList.end(); ++i)
    {
      aVoxel->surfaceVoxels->push_back(*i);
    }
  aVoxel->adjoiningSize = forward-aVoxel->adjoiningVoxels;
}

Species* SpatiocyteStepper::id2species(unsigned short id)
{
  return theSpecies[id];
}

Comp* SpatiocyteStepper::id2Comp(unsigned short id)
{
  return theSpecies[id]->getComp();
}

Voxel* SpatiocyteStepper::coord2voxel(unsigned int aCoord)
{
  return &theLattice[aCoord];
}

Comp* SpatiocyteStepper::system2Comp(System* aSystem)
{
  for(unsigned int i(0); i != theComps.size(); ++i)
    {
      if(theComps[i]->system == aSystem)
        {
          return theComps[i];
        }
    }
  return NULL;
}

void SpatiocyteStepper::setSurfaceSubunit(Voxel* aVoxel,
                                          Comp* aComp)
{
  // The subunit is only useful for a cylindrical surface
  // and for polymerization on it.
  aVoxel->subunit = new Subunit;
  aVoxel->subunit->voxel = aVoxel;
  Point& aPoint(aVoxel->subunit->surfacePoint);
  aPoint = coord2point(aVoxel->coord);
  double aRadius(aComp->lengthY/2);
  Point aCenterPoint(aComp->centerPoint);
  Point aWestPoint(aComp->centerPoint);
  Point anEastPoint(aComp->centerPoint); 
  aWestPoint.x = aComp->centerPoint.x-aComp->lengthX/2+aComp->lengthY/2;
  anEastPoint.x = aComp->centerPoint.x+aComp->lengthX/2-aComp->lengthY/2;
  if(aPoint.x < aWestPoint.x)
    {
      aCenterPoint.x = aWestPoint.x;
    }
  else if(aPoint.x > anEastPoint.x)
    {
      aCenterPoint.x = anEastPoint.x;
    }
  else
    {
      aCenterPoint.x = aPoint.x;
    }
  double X(aPoint.x-aCenterPoint.x);
  double Y(aPoint.y-aCenterPoint.y);
  double Z(aPoint.z-aCenterPoint.z);
  double f(atan2(X, Z));
  double d(atan2(sqrt(X*X+Z*Z),Y));
  aPoint.x = aCenterPoint.x + sin(f)*aRadius*sin(d);
  aPoint.y = aCenterPoint.y + aRadius*cos(d);
  aPoint.z = aCenterPoint.z + cos(f)*aRadius*sin(d);
}

void SpatiocyteStepper::setIntersectingCompartmentList() 
{
  setIntersectingPeers();
  setIntersectingParent();
}

void SpatiocyteStepper::setIntersectingPeers()
{
  for(std::vector<Comp*>::reverse_iterator i(theComps.rbegin());
      i != theComps.rend(); ++i)
    {
      //Only proceed if the volume is not enclosed:
      if(!(*i)->system->isRootSystem() && (*i)->dimension == 3)
        {
          for(std::vector<Comp*>::reverse_iterator j(theComps.rbegin());
              j != theComps.rend(); ++j)
            {
              //Only proceed if the volume of the peer is enclosed:
              if((*i) != (*j) && !(*j)->system->isRootSystem() &&
                 (*j)->dimension == 3)
                {
                  //If i and j are peer volume compartments:
                  if((*i)->system->getSuperSystem() == 
                     (*j)->system->getSuperSystem())
                    {
                      Point a((*i)->centerPoint);
                      Point b((*j)->centerPoint);
                      if(((a.x+(*i)->lengthX/2 > b.x-(*j)->lengthX/2) ||
                          (a.x-(*i)->lengthX/2 < b.x+(*j)->lengthX/2)) && 
                         ((a.y+(*i)->lengthY/2 > b.y-(*j)->lengthY/2) ||
                          (a.y-(*i)->lengthY/2 < b.y+(*j)->lengthY/2)) &&
                         ((a.z+(*i)->lengthZ/2 > b.z-(*j)->lengthZ/2) ||
                          (a.z-(*i)->lengthZ/2 < b.z+(*j)->lengthZ/2)))
                        {
                          if((*i)->enclosed <= (*j)->enclosed) 
                            {
                              (*i)->intersectPeers.push_back(*j);
                            }
                          else
                            {
                              (*i)->intersectLowerPeers.push_back(*j);
                            }
                        }
                    }
                }
            }
        }
    }
}

void SpatiocyteStepper::setIntersectingParent() 
{
  for(std::vector<Comp*>::iterator i(theComps.begin());
      i != theComps.end(); ++i)
    {
      (*i)->isIntersectRoot = false;
      (*i)->isIntersectParent = false;
      //only proceed if the volume is not enclosed:
      if(!(*i)->system->isRootSystem() && (*i)->dimension == 3)
        {
          Comp* aParentComp(system2Comp((*i)->system->getSuperSystem())); 
          Point a((*i)->centerPoint);
          Point b(aParentComp->centerPoint);
          if((a.x+(*i)->lengthX/2 > b.x+aParentComp->lengthX/2) ||
             (a.x-(*i)->lengthX/2 < b.x-aParentComp->lengthX/2) || 
             (a.y+(*i)->lengthY/2 > b.y+aParentComp->lengthY/2) ||
             (a.y-(*i)->lengthY/2 < b.y-aParentComp->lengthY/2) ||
             (a.z+(*i)->lengthZ/2 > b.z+aParentComp->lengthZ/2) ||
             (a.z-(*i)->lengthZ/2 < b.z-aParentComp->lengthZ/2))
            {
              if(aParentComp->system->isRootSystem())
                {
                  (*i)->isIntersectRoot = true;
                }
              else
                {
                  (*i)->isIntersectParent = true;
                }
            }
        }
    }
}

void SpatiocyteStepper::compartmentalizeLattice() 
{
  for(std::vector<Voxel>::iterator i(theLattice.begin());
      i != theLattice.end(); ++i)
    {
      if((*i).id != theNullID)
        { 
          compartmentalizeVoxel(&(*i), theComps[0]);
        }
    }
}

bool SpatiocyteStepper::compartmentalizeVoxel(Voxel* aVoxel, Comp* aComp)
{
  if(aComp->dimension == 3)
    {
      if(aComp->system->isRootSystem() ||
         isInsideCoord(aVoxel->coord, aComp, 0))
        {
          //Check if the voxel belongs to an intersecting peer compartment:
          if(!aComp->intersectPeers.empty())
            {
              //The function isPeerVoxel also checks if the voxel is a future
              //surface voxel of the peer compartment: 
              if(isPeerVoxel(aVoxel, aComp))
                { 
                  return false;
                }
              if(aComp->surfaceSub && aComp->surfaceSub->enclosed)
                {
                  //Check if the voxel is neighbor of a peer voxel (either
                  //a future surface voxel)
                  if(isEnclosedSurfaceVoxel(aVoxel, aComp))
                    {
                      aVoxel->id = aComp->surfaceSub->vacantID;
                      aComp->surfaceSub->coords.push_back(
                                                 aVoxel->coord-theStartCoord);
                      setMinMaxSurfaceDimensions(aVoxel->coord, aComp);
                      return true;
                    }
                }
            }
          if(aComp->isIntersectRoot)
            {
              Comp* aRootComp(system2Comp(aComp->system->getSuperSystem())); 
              if(aRootComp->surfaceSub && 
                 aRootComp->surfaceSub->enclosed && 
                 isRootSurfaceVoxel(aVoxel, aRootComp))
                {
                  return false;
                }
              if(aComp->surfaceSub && 
                 isEnclosedRootSurfaceVoxel(aVoxel, aComp, aRootComp))
                {
                  aVoxel->id = aComp->surfaceSub->vacantID;
                  aComp->surfaceSub->coords.push_back(
                                                aVoxel->coord-theStartCoord);
                  setMinMaxSurfaceDimensions(aVoxel->coord, aComp);
                  return true;
                }
            }
          else if(aComp->isIntersectParent)
            {
              Comp* aParentComp(system2Comp(aComp->system->getSuperSystem())); 
              if(aComp->surfaceSub && aComp->surfaceSub->enclosed &&
                 isParentSurfaceVoxel(aVoxel, aParentComp))
                {
                  aVoxel->id = aComp->surfaceSub->vacantID;
                  aComp->surfaceSub->coords.push_back(
                                                aVoxel->coord-theStartCoord);
                  setMinMaxSurfaceDimensions(aVoxel->coord, aComp);
                  return true;
                }
            }
          for(unsigned int i(0); i != aComp->immediateSubs.size(); ++i)
            {
              if(compartmentalizeVoxel(aVoxel, aComp->immediateSubs[i]))
                {
                  return true;
                }
            }
          aVoxel->id = aComp->vacantID;
          aComp->coords.push_back(aVoxel->coord-theStartCoord);
          return true;
        }
      if(aComp->surfaceSub)
        { 
          if(isInsideCoord(aVoxel->coord, aComp, 4) &&
             isSurfaceVoxel(aVoxel, aComp))
            {
              aVoxel->id = aComp->surfaceSub->vacantID;
              aComp->surfaceSub->coords.push_back(
                                              aVoxel->coord-theStartCoord);
              setMinMaxSurfaceDimensions(aVoxel->coord, aComp);
              return true;
            }
        }
    }
  else if(aComp->system->getSuperSystem()->isRootSystem())
    {
      Comp* aRootComp(system2Comp(aComp->system->getSuperSystem())); 
      if(!isInsideCoord(aVoxel->coord, aRootComp, -4) &&
         isRootSurfaceVoxel(aVoxel, aRootComp))
        {
          aVoxel->id = aComp->vacantID;
          aComp->coords.push_back(aVoxel->coord-theStartCoord);
          setMinMaxSurfaceDimensions(aVoxel->coord, aRootComp);
          return true;
        }
    }
  return false;
}

bool SpatiocyteStepper::isRootSurfaceVoxel(Voxel* aVoxel, Comp* aComp)
{
  for(unsigned int i(0); i != theAdjoiningVoxelSize; ++i)
    {
      if(aVoxel->adjoiningVoxels[i]->id == theNullID ||
         aVoxel->adjoiningVoxels[i] == aVoxel)
        {
          return true;
        }
    }
  return false;
}

bool SpatiocyteStepper::isParentSurfaceVoxel(Voxel* aVoxel, Comp* aComp)
{
  if(!isInsideCoord(aVoxel->coord, aComp, -4))
    {
      for(unsigned int i(0); i != theAdjoiningVoxelSize; ++i)
        {
          if(!isInsideCoord(aVoxel->adjoiningVoxels[i]->coord, aComp, 0))
            {
              return true;
            }
        }
    }
  return false;
}

bool SpatiocyteStepper::isEnclosedRootSurfaceVoxel(Voxel* aVoxel, Comp* aComp,
                                                   Comp* aRootComp)
{ 
  if(aComp->surfaceSub->enclosed || aRootComp->enclosed <= aComp->enclosed)
    { 
      if(!isInsideCoord(aVoxel->coord, aRootComp, -4))
        {
          if(!aRootComp->surfaceSub || 
             (aRootComp->surfaceSub && !aRootComp->surfaceSub->enclosed))
            {
              if(isRootSurfaceVoxel(aVoxel, aComp))
                {
                  return true;
                }
            }
          else
            {
              for(unsigned int i(0); i != theAdjoiningVoxelSize; ++i)
                {
                  if(isRootSurfaceVoxel(aVoxel->adjoiningVoxels[i], aComp))
                    {
                      return true;
                    }
                }
            }
        }
    }
  return false;
}

bool SpatiocyteStepper::isSurfaceVoxel(Voxel* aVoxel, Comp* aComp)
{
  if(isPeerVoxel(aVoxel, aComp))
    { 
      return false;
    }
  if(aComp->surfaceSub && !aComp->surfaceSub->enclosed &&
     isLowerPeerVoxel(aVoxel, aComp))
    {
      return false;
    }
  if(aComp->isIntersectRoot)
    {
      Comp* aRootComp(system2Comp(aComp->system->getSuperSystem())); 
      if(aRootComp->surfaceSub && aRootComp->surfaceSub->enclosed && 
         isRootSurfaceVoxel(aVoxel, aRootComp))
        {
          return false;
        }
    }
  for(unsigned int i(0); i != theAdjoiningVoxelSize; ++i)
    {
      if(isInsideCoord(aVoxel->adjoiningVoxels[i]->coord, aComp, 0))
        {
          return true;
        }
    }
  return false;
}

bool SpatiocyteStepper::isEnclosedSurfaceVoxel(Voxel* aVoxel, Comp* aComp)
{
  for(unsigned int i(0); i != theAdjoiningVoxelSize; ++i)
    {
      if(isPeerVoxel(aVoxel->adjoiningVoxels[i], aComp))
        {
          return true;
        }
    }
  return false;
}

bool SpatiocyteStepper::isLowerPeerVoxel(Voxel* aVoxel, Comp* aComp)
{
  for(std::vector<Comp*>::iterator i(aComp->intersectLowerPeers.begin());
      i != aComp->intersectLowerPeers.end(); ++i)
    {
      if(isInsideCoord(aVoxel->coord, *i, 0))
        {
          return true;
        }
    }
  return false;
}

bool SpatiocyteStepper::isPeerVoxel(Voxel* aVoxel, Comp* aComp)
{
  for(std::vector<Comp*>::iterator i(aComp->intersectPeers.begin());
      i != aComp->intersectPeers.end(); ++i)
    {
      if(isInsideCoord(aVoxel->coord, *i, 0))
        {
          return true;
        }
    }
  //The voxel is not inside one of the intersecting peer compartments,
  //we are now checking if the one of voxel's neighbor belongs to the
  //peer comp, to determine if aVoxel is a surface voxel of the peer:
  for(std::vector<Comp*>::iterator i(aComp->intersectPeers.begin());
      i != aComp->intersectPeers.end(); ++i)
    {
      if((*i)->surfaceSub && (*i)->surfaceSub->enclosed)
        { 
          if(isInsideCoord(aVoxel->coord, *i, 4))
            {
              for(unsigned int j(0); j != theAdjoiningVoxelSize; ++j)
                {
                  if(isInsideCoord(aVoxel->adjoiningVoxels[j]->coord, *i, 0))
                    {
                      return true;
                    }
                }
            }
        }
    }
  return false;
}

bool SpatiocyteStepper::isLineVoxel(Voxel* aVoxel, Comp* aComp)
{
    const double safety(2.0);
    const Point aPoint(coord2point(aVoxel->coord));

    double distance(aPoint.x - aComp->centerPoint.x);
    if (-safety < distance && distance <= 0)
    {
        // This is not efficient because we don't need to check volume voxels.
        // However, at this time, aVoxel->adjoinigVoxels is not properly 
        // aligned yet.
        for(unsigned int i(0); i != theAdjoiningVoxelSize; ++i)
        {
            const Voxel* adjoiningVoxel(aVoxel->adjoiningVoxels[i]);
            const Point adjoiningPoint(coord2point(adjoiningVoxel->coord));
            const double distance_i(adjoiningPoint.x - aComp->centerPoint.x);
            if (distance_i > 0)
            {
                return true;
            }
        }
    }
                       
    distance = aPoint.y - aComp->centerPoint.y;
    if (-safety < distance && distance <= 0)
    {
        // This is not efficient because we don't need to check volume voxels.
        // However, at this time, aVoxel->adjoinigVoxels is not properly 
        // aligned yet.
        for(unsigned int i(0); i != theAdjoiningVoxelSize; ++i)
        {
            const Voxel* adjoiningVoxel(aVoxel->adjoiningVoxels[i]);
            const Point adjoiningPoint(coord2point(adjoiningVoxel->coord));
            const double distance_i(adjoiningPoint.y - aComp->centerPoint.y);
            if (distance_i > 0)
            {
                return true;
            }
        }
    }
    return false;
}



void SpatiocyteStepper::setMinMaxSurfaceDimensions(unsigned int aCoord, 
                                                   Comp* aComp)
{
  unsigned int aRow;
  unsigned int aLayer;
  unsigned int aCol;
  coord2global(aCoord, &aRow, &aLayer, &aCol);
  if(aRow < aComp->minRow)
    {
      aComp->minRow = aRow;
      aComp->surfaceSub->minRow = aRow;
    }
  else if(aRow > aComp->maxRow)
    {
      aComp->maxRow = aRow;
      aComp->surfaceSub->maxRow = aRow;
    }
  if(aCol < aComp->minCol)
    {
      aComp->minCol = aCol;
      aComp->surfaceSub->minCol = aCol;
    }
  else if(aCol > aComp->maxCol)
    {
      aComp->maxCol = aCol;
      aComp->surfaceSub->maxCol = aCol;
    }
  if(aLayer < aComp->minLayer)
    {
      aComp->minLayer = aLayer;
      aComp->surfaceSub->minLayer = aLayer;
    }
  else if(aLayer > aComp->maxLayer)
    {
      aComp->maxLayer = aLayer;
      aComp->surfaceSub->maxLayer = aLayer;
    }
}



void SpatiocyteStepper::rotateX(double angle, Point* aPoint)
{ 
  if(angle)
    {
      double y(aPoint->y);
      double z(aPoint->z);
      aPoint->y = y*cos(angle)-z*sin(angle);
      aPoint->z = y*sin(angle)+z*cos(angle);
    }
}

void SpatiocyteStepper::rotateY(double angle, Point* aPoint)
{ 
  if(angle)
    {
      double x(aPoint->x);
      double z(aPoint->z);
      aPoint->x = x*cos(angle)+z*sin(angle);
      aPoint->z = z*cos(angle)-x*sin(angle);
    }
}

void SpatiocyteStepper::rotateZ(double angle, Point* aPoint)
{ 
  if(angle)
    {
      double x(aPoint->x);
      double y(aPoint->y);
      aPoint->x = x*cos(angle)-y*sin(angle);
      aPoint->y = x*sin(angle)+y*cos(angle);
    }
}

bool SpatiocyteStepper::isInsideCoord(unsigned int aCoord,
                                      Comp* aComp, double delta)
{
  Point aPoint(coord2point(aCoord));
  Point aCenterPoint(aComp->centerPoint);
  Point aWestPoint(aComp->centerPoint);
  Point anEastPoint(aComp->centerPoint); 
  aPoint.x = aPoint.x - aCenterPoint.x;
  aPoint.y = aPoint.y - aCenterPoint.y;
  aPoint.z = aPoint.z - aCenterPoint.z;
  rotateX(aComp->rotateX, &aPoint);
  rotateY(aComp->rotateY, &aPoint);
  rotateZ(aComp->rotateZ, &aPoint);
  aPoint.x = aPoint.x + aCenterPoint.x;
  aPoint.y = aPoint.y + aCenterPoint.y;
  aPoint.z = aPoint.z + aCenterPoint.z;
  double aRadius(0);
  switch(aComp->geometry)
    {
    case CUBOID:
      if(sqrt(pow(aPoint.x-aCenterPoint.x, 2)) <= 
         aComp->lengthX/2+theNormalizedVoxelRadius+delta &&
         sqrt(pow(aPoint.y-aCenterPoint.y, 2)) <= 
         aComp->lengthY/2+theNormalizedVoxelRadius+delta &&
         sqrt(pow(aPoint.z-aCenterPoint.z, 2)) <= 
         aComp->lengthZ/2+theNormalizedVoxelRadius+delta)
        {
          return true;
        }
      break;
    case ELLIPSOID:
      //If the distance between the voxel and the center point is less than 
      //or equal to radius-2, then the voxel cannot be a surface voxel:
      if(pow(aPoint.x-aCenterPoint.x, 2)/pow((aComp->lengthX+delta)/2, 2)+ 
         pow(aPoint.y-aCenterPoint.y, 2)/pow((aComp->lengthY+delta)/2, 2)+ 
         pow(aPoint.z-aCenterPoint.z, 2)/pow((aComp->lengthZ+delta)/2, 2) <= 1)
        {
          return true;
        }
      break;
    case CYLINDER: 
      //The axial point of the cylindrical portion of the rod:
      aCenterPoint.x = aPoint.x; 
      aWestPoint.x = aComp->centerPoint.x-aComp->lengthX/2;
      anEastPoint.x = aComp->centerPoint.x+aComp->lengthX/2;
      aRadius = aComp->lengthY/2+theNormalizedVoxelRadius;
      //If the distance between the voxel and the center point is less than 
      //or equal to the radius, then the voxel must be inside the Comp:
      if((aPoint.x >= aWestPoint.x && aPoint.x <= anEastPoint.x &&
          getDistance(&aPoint, &aCenterPoint) <= aRadius+delta))
        { 
          return true;
        }
      break;
    case ROD: 
      //The axial point of the cylindrical portion of the rod:
      aCenterPoint.x = aPoint.x; 
      aWestPoint.x = aComp->centerPoint.x-aComp->lengthX/2+aComp->lengthY/2;
      anEastPoint.x = aComp->centerPoint.x+aComp->lengthX/2-aComp->lengthY/2;
      aRadius = aComp->lengthY/2+theNormalizedVoxelRadius;
      //If the distance between the voxel and the center point is less than 
      //or equal to the radius, then the voxel must be inside the Comp:
      if((aPoint.x >= aWestPoint.x && aPoint.x <= anEastPoint.x &&
          getDistance(&aPoint, &aCenterPoint) <= aRadius+delta) ||
         (aPoint.x < aWestPoint.x &&
          getDistance(&aPoint, &aWestPoint) <= aRadius+delta) ||
         (aPoint.x > anEastPoint.x &&
          getDistance(&aPoint, &anEastPoint) <= aRadius+delta))
        { 
          return true;
        }
      break;
    case TORUS: 
      return true;
      break;
    case PYRAMID: 
      aRadius = ((aCenterPoint.y+aComp->lengthY/2)-aPoint.y)/aComp->lengthY;
      if(sqrt(pow(aPoint.y-aCenterPoint.y, 2)) <= aComp->lengthY/2+delta &&
         sqrt(pow(aPoint.x-aCenterPoint.x, 2)) <= 
         aComp->lengthX*aRadius/2+delta &&
         sqrt(pow(aPoint.z-aCenterPoint.z, 2)) <=
         aComp->lengthZ*aRadius/2+delta)
        {
          return true;
        }
      break;
    }
  return false;
}

void SpatiocyteStepper::populateComp(Comp* aComp)
{
  unsigned int populationSize(0);
  std::vector<Species*> diffuseVacantSpecies;
  std::vector<Species*> diffuseSpecies;
  std::vector<Species*> normalSpecies;
  for(std::vector<Species*>::const_iterator i(aComp->species.begin());
      i != aComp->species.end(); ++i)
    {
      if((*i)->getIsDiffuseVacant())
        {
          populationSize += (unsigned int)(*i)->getPopulateMoleculeSize();
          diffuseVacantSpecies.push_back(*i);
        }
      else if((*i)->getVacantSpecies()->getIsDiffuseVacant())
        {
          diffuseSpecies.push_back(*i);
        }
      else
        {
          populationSize += (unsigned int)(*i)->getPopulateMoleculeSize();
          normalSpecies.push_back(*i);
        }
    }
  if(double(populationSize)/aComp->coords.size() > 0.2)
    { 
      populateSpeciesDense(diffuseVacantSpecies, populationSize,
                           aComp->coords.size());
      populateSpeciesDense(normalSpecies, populationSize, aComp->coords.size());
    }
  else
    {
      populateSpeciesSparse(diffuseVacantSpecies);
      populateSpeciesSparse(normalSpecies);
    }
  for(std::vector<Species*>::const_iterator i(diffuseSpecies.begin());
      i != diffuseSpecies.end(); ++i)
    {
      (*i)->populateUniformDiffuseVacant();
    }
}

void SpatiocyteStepper::populateSpeciesDense(std::vector<Species*>&
                                             aSpeciesList, unsigned int aSize,
                                             unsigned int availableVoxelSize)
{
  unsigned int count(0);
  unsigned int* populateVoxels(new unsigned int[aSize]);
  unsigned int* availableVoxels(new unsigned int [availableVoxelSize]); 
  for(unsigned int i(0); i != availableVoxelSize; ++i)
    {
      availableVoxels[i] = i;
    }
  gsl_ran_choose(getRng(), populateVoxels, aSize, availableVoxels,
                 availableVoxelSize, sizeof(unsigned int));
  //gsl_ran_choose arranges the position ascending, so we need
  //to shuffle the order of voxel positions:
  gsl_ran_shuffle(getRng(), populateVoxels, aSize, sizeof(unsigned int)); 
  for(std::vector<Species*>::const_iterator i(aSpeciesList.begin());
      i != aSpeciesList.end(); ++i)
    {
      (*i)->populateCompUniform(populateVoxels, &count);
    }
  delete[] populateVoxels;
  delete[] availableVoxels;
}

void SpatiocyteStepper::populateSpeciesSparse(std::vector<Species*>&
                                              aSpeciesList)
{
  for(std::vector<Species*>::const_iterator i(aSpeciesList.begin());
      i != aSpeciesList.end(); ++i)
    {
      (*i)->populateCompUniformSparse();
    }
}


void SpatiocyteStepper::clearComp(Comp* aComp)
{
  for(std::vector<Species*>::const_iterator i(aComp->species.begin());
      i != aComp->species.end(); ++i)
    {
      (*i)->removeMolecules();
    }
}


std::vector<Comp*> const& SpatiocyteStepper::getComps() const
{
  return theComps;
}

