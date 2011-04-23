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
  //We need a compartment tree to assign the voxels to each compartment
  //and get the available number of vacant voxels. The compartmentalized
  //vacant voxels are needed to randomly place molecules according to the
  //compartment:
  std::cout << "2. creating compartments..." << std::endl;
  registerCompartments();
  setCompartmentsProperties();
  std::cout << "3. setting up lattice properties..." << std::endl;
  setLatticeProperties(); 
  setCompartmentsCenterPoint();
  //All species have been created at this point, we initialize them now:
  std::cout << "4. initializing species..." << std::endl;
  initSpecies();
  std::cout << "5. initializing processes the second time..." << std::endl;
  initProcessSecond();
  std::cout << "7. constructing lattice..." << std::endl;
  constructLattice();
  std::cout << "8. shuffling adjoining voxels..." << std::endl;
  shuffleAdjoiningVoxels();
  std::cout << "9. compartmentalizing lattice..." << std::endl;
  compartmentalizeLattice();
  std::cout << "10. setting up properties of surface voxels..." << std::endl;
  setSurfaceVoxelProperties();
  std::cout << "11. populating compartments with molecules..." << std::endl;
  populateCompartments();
  storeSimulationParameters();
  //checkSurfaceCompartment();
  std::cout << "12. initializing processes the third time..." << std::endl;
  initProcessThird();
  std::cout << "13. initializing the priority queue..." << std::endl;
  initPriorityQueue();
  printSimulationParameters();
  std::cout << "14. initializing processes the fourth time..." << std::endl;
  initProcessFourth();
  initProcessLastOnce();
  std::cout << "15. simulation is started..." << std::endl;
  //checkSurfaceCompartment();
  //checkLattice();
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
  clearCompartments();
  populateCompartments();
  initProcessThird();
  initPriorityQueue();
  initProcessFourth();
  //checkLattice();
}

Species* SpatiocyteStepper::addSpecies(Variable* aVariable)
{
  std::vector<Species*>::iterator aSpeciesIter(variable2species(aVariable));
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
  std::vector<Species*>::iterator aSpeciesIter(variable2species(aVariable));
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
      //If the voxel is on one of the 6 cubic surfaces:
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
SpatiocyteStepper::variable2species(Variable* aVariable)
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

void SpatiocyteStepper::checkSurfaceCompartment()
{
  for(unsigned int i(0); i!=theLattice.size(); ++i)
    {
      Voxel* aVoxel(&theLattice[i]);
      for(int k(0); k != ADJOINING_VOXEL_SIZE; ++k)
        {
          if(aVoxel == aVoxel->adjoiningVoxels[k])
            {
              theSpecies[1]->addSimpleMolecule(aVoxel);
              break;
            }
        }
    }
}

/*
void SpatiocyteStepper::checkSurfaceCompartment()
{
  for(vector<Compartment*>::const_iterator i(theCompartments.begin());
      i != theCompartments.end(); ++i)
    {
      if((*i)->isSurface)
        {
          for(vector<unsigned int>::const_iterator j((*i)->coords.begin());
              j != (*i)->coords.end(); ++j )
            {
              Voxel* aVoxel(&theLattice[*j]);
              for(int k(0); k != ADJOINING_VOXEL_SIZE; ++k)
                {
                  if(aVoxel == aVoxel->adjoiningVoxels[k])
                    {
                      theSpecies[2]->addSimpleMolecule(aVoxel);
                      break;
                    }
                }
            }
        }
    }
}
*/

/*
void SpatiocyteStepper::checkSurfaceCompartment()
{
  vector<int> surfaceCnt;
  surfaceCnt.resize(12);
  for(int i(0); i!=12; ++i)
    {
      surfaceCnt[i] = 0;
    }
  for(vector<Compartment*>::const_iterator i(theCompartments.begin());
      i != theCompartments.end(); ++i)
    {
      if((*i)->isSurface)
        {
          std::cout << "size:" << (*i)->coords.size() << std::endl;
          int surfaceID((*i)->vacantID);
          for(vector<unsigned int>::const_iterator j((*i)->coords.begin());
              j != (*i)->coords.end(); ++j )
            {
              Voxel* aVoxel(&theLattice[*j]);
              int cnt(0);
              for(int k(0); k != ADJOINING_VOXEL_SIZE; ++k)
                {
                  if(aVoxel != aVoxel->adjoiningVoxels[k] &&
                     aVoxel->adjoiningVoxels[k]->id == surfaceID)
                    {
                      ++cnt;
                    }
                }
              ++surfaceCnt[cnt];
              if(cnt > 2 && cnt < 9)
                {
                  theSpecies[theSpecies.size()-9+cnt]->addSimpleMolecule(aVoxel);
                }
            }
        }
    }
  int total(0);
  for(int i(0); i!=12; ++i)
    {
      total += surfaceCnt[i];
      std::cout << i << ": " << surfaceCnt[i] << std::endl;
    }
  std::cout << "total:" << total << std::endl;
  std::cout << "size:" << theSpecies.back()->size() << std::endl;
}
*/

void SpatiocyteStepper::initSpecies()
{
  for(std::vector<Species*>::iterator i(theSpecies.begin());
      i != theSpecies.end(); ++i)
    {
      (*i)->initialize(theSpecies.size());
    }
}

void SpatiocyteStepper::initProcessSecond()
{
  for(std::vector<Process*>::const_iterator i(theProcessVector.begin());
      i != theProcessVector.end(); ++i)
    {      
      SpatiocyteProcessInterface* aProcess(dynamic_cast<SpatiocyteProcessInterface*>(*i));
      aProcess->initializeSecond();
    }
}

void SpatiocyteStepper::printProcessParameters()
{
  for(std::vector<Process*>::const_iterator i(theProcessVector.begin());
      i != theProcessVector.end(); ++i)
    {      
      SpatiocyteProcessInterface* aProcess(dynamic_cast<SpatiocyteProcessInterface*>(*i));
      aProcess->printParameters();
    }
}

void SpatiocyteStepper::initProcessThird()
{
  for(std::vector<Process*>::const_iterator i(theProcessVector.begin());
      i != theProcessVector.end(); ++i)
    {      
      SpatiocyteProcessInterface* aProcess(dynamic_cast<SpatiocyteProcessInterface*>(*i));
      aProcess->initializeThird();
    }
}

void SpatiocyteStepper::initProcessFourth()
{
  for(std::vector<Process*>::const_iterator i(theProcessVector.begin());
      i != theProcessVector.end(); ++i)
    {      
      SpatiocyteProcessInterface* aProcess(dynamic_cast<SpatiocyteProcessInterface*>(*i));
      aProcess->initializeFourth();
    }
  setStepInterval(thePriorityQueue.getTop()->getTime()-getCurrentTime());
}

void SpatiocyteStepper::initProcessLastOnce()
{
  for(std::vector<Process*>::const_iterator i(theProcessVector.begin());
      i != theProcessVector.end(); ++i)
    {      
      SpatiocyteProcessInterface* aProcess(dynamic_cast<SpatiocyteProcessInterface*>(*i));
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
             aClassName == "MoleculePopulateProcess" ||
             aClassName == "CoordinateLogProcess" ||
             aClassName == "VisualizationLogProcess" ||
             aClassName == "FluorescentImagingProcess" ||
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

void SpatiocyteStepper::populateCompartments()
{
  for(std::vector<Compartment*>::const_iterator i(theCompartments.begin());
      i != theCompartments.end(); ++i)
    {
      populateCompartment(*i);
    }
}

void SpatiocyteStepper::clearCompartments()
{
  for(std::vector<Compartment*>::const_iterator i(theCompartments.begin());
      i != theCompartments.end(); ++i)
    {
      clearCompartment(*i);
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


void SpatiocyteStepper::registerCompartments()
{
  System* aRootSystem(getModel()->getRootSystem());
  std::vector<Compartment*> allSubs;
  //The root compartment is theCompartments[0]
  theCompartments.push_back(registerCompartment(aRootSystem, &allSubs));
  //After this we will create an species to get an ID to represent
  //NULL compartments. So let us specify the correct
  //size of the biochemical species to be simulated before additional
  //non-biochemical species are created:
  theBioSpeciesSize = theSpecies.size();
  //Create one last species to represent a NULL compartment. This is for
  //voxels that do not belong to any compartments:
  Species* aSpecies(new Species(this, NULL, theSpecies.size(), 0, getRng()));
  theSpecies.push_back(aSpecies);
  aSpecies->setCompartment(NULL);
  theNullID = aSpecies->getID(); 
  //Expand the tree of immediate subcompartments into single list such that
  //the super compartments come first while the subcompartments 
  //come later in the list:
  std::vector<Compartment*> compartments(theCompartments[0]->immediateSubs);
  while(!compartments.empty())
    {
      std::vector<Compartment*> subCompartments;
      for(unsigned int i(0); i != compartments.size(); ++i)
        {
          theCompartments.push_back(compartments[i]);
          for(unsigned int j(0);
              j != compartments[i]->immediateSubs.size(); ++j)
            {
              subCompartments.push_back(compartments[i]->immediateSubs[j]);
            }
        }
      compartments = subCompartments;
    }
}

//allSubs contains all the subCompartments (child, grand child, great grand
//child, etc). Used to calculate the total number of compartment voxels.
Compartment* SpatiocyteStepper::registerCompartment(System* aSystem,
                                            std::vector<Compartment*>* allSubs)
{ 
  //We execute this function to register the System, and its subsystems
  //recursively.
  Compartment* aCompartment(new Compartment);
  aCompartment->lengthX = 0;
  aCompartment->lengthY = 0;
  aCompartment->lengthZ = 0;
  aCompartment->originX = 0;
  aCompartment->originY = 0;
  aCompartment->originZ = 0;
  aCompartment->rotateX = 0;
  aCompartment->rotateY = 0;
  aCompartment->rotateZ = 0;
  aCompartment->xyPlane = 0;
  aCompartment->xzPlane = 0;
  aCompartment->yzPlane = 0;
  aCompartment->specVolume = 0;
  aCompartment->system = aSystem;
  aCompartment->surfaceSub = NULL;
  //Default compartment shape is spherical:
  aCompartment->shape = 0;
  //Default is volume compartment:
  aCompartment->isSurface = false;
  if(getVariable(aSystem, "TYPE"))
    { 
      if(aSystem->getVariable("TYPE")->getValue() == SURFACE)
        {
          aCompartment->isSurface = true;
        }
    }
  if(!aCompartment->isSurface)
    {
      if(getVariable(aSystem, "SHAPE"))
        { 
          aCompartment->shape = aSystem->getVariable("SHAPE")->getValue();
        }
      if(getVariable(aSystem, "LENGTHX"))
        {
          aCompartment->lengthX = aSystem->getVariable("LENGTHX")->getValue();
        }
      if(getVariable(aSystem, "LENGTHY"))
        {
          aCompartment->lengthY = aSystem->getVariable("LENGTHY")->getValue();
        }
      if(getVariable(aSystem, "LENGTHZ"))
        {
          aCompartment->lengthZ = aSystem->getVariable("LENGTHZ")->getValue();
        }
      if(getVariable(aSystem, "ORIGINX"))
        {
          aCompartment->originX = aSystem->getVariable("ORIGINX")->getValue();
        }
      if(getVariable(aSystem, "ORIGINY"))
        {
          aCompartment->originY = aSystem->getVariable("ORIGINY")->getValue();
        }
      if(getVariable(aSystem, "ORIGINZ"))
        {
          aCompartment->originZ = aSystem->getVariable("ORIGINZ")->getValue();
        }
      if(getVariable(aSystem, "ROTATEX"))
        {
          aCompartment->rotateX = aSystem->getVariable("ROTATEX")->getValue();
        }
      if(getVariable(aSystem, "ROTATEY"))
        {
          aCompartment->rotateY = aSystem->getVariable("ROTATEY")->getValue();
        }
      if(getVariable(aSystem, "ROTATEZ"))
        {
          aCompartment->rotateZ = aSystem->getVariable("ROTATEZ")->getValue();
        }
      if(getVariable(aSystem, "XYPLANE"))
        {
          aCompartment->xyPlane = aSystem->getVariable("XYPLANE")->getValue();
        }
      if(getVariable(aSystem, "XZPLANE"))
        {
          aCompartment->xzPlane = aSystem->getVariable("XZPLANE")->getValue();
        }
      if(getVariable(aSystem, "YZPLANE"))
        {
          aCompartment->yzPlane = aSystem->getVariable("YZPLANE")->getValue();
        }
      if(getVariable(aSystem, "SIZE"))
        { 
          aCompartment->specVolume = aSystem->getVariable("SIZE")->getValue();
          //Change SIZE unit to liter to be consistent with E-Cell's SIZE unit.
          aSystem->getVariable("SIZE")->setValue(aCompartment->specVolume*1e+3);
        }
    }
  registerCompartmentSpecies(aCompartment);
  //Systems contains all the subsystems of a System.
  //For example /membrane is the subsystem of /:
  FOR_ALL(System::Systems, aSystem->getSystems())
    {
      Compartment* aSubCompartment(registerCompartment(i->second, allSubs)); 
      allSubs->push_back(aSubCompartment);
      aCompartment->immediateSubs.push_back(aSubCompartment);
      if(aSubCompartment->isSurface)
        {
          aSubCompartment->shape = aCompartment->shape;
          aSubCompartment->specVolume = aCompartment->specVolume;
          aSubCompartment->lengthX = aCompartment->lengthX;
          aSubCompartment->lengthY = aCompartment->lengthY;
          aSubCompartment->lengthZ = aCompartment->lengthZ;
          aSubCompartment->originX = aCompartment->originX;
          aSubCompartment->originY = aCompartment->originY;
          aSubCompartment->originZ = aCompartment->originZ;
          aSubCompartment->xyPlane = aCompartment->xyPlane;
          aSubCompartment->xzPlane = aCompartment->xzPlane;
          aSubCompartment->yzPlane = aCompartment->yzPlane;
          aCompartment->surfaceSub = aSubCompartment;
        }
    }
  aCompartment->allSubs = *allSubs;
  return aCompartment;
}

void SpatiocyteStepper::setCompartmentsProperties()
{
  for(unsigned int i(0); i != theCompartments.size(); ++i)
    {
      std::cout << theCompartments[i]->system->getFullID().asString() << std::endl;
      setCompartmentProperties(theCompartments[i]);
    }
}

void SpatiocyteStepper::setCompartmentsCenterPoint()
{
  for(unsigned int i(0); i != theCompartments.size(); ++i)
    {
      setCompartmentCenterPoint(theCompartments[i]);
    }
}

void SpatiocyteStepper::registerCompartmentSpecies(Compartment* aCompartment)
{
  System* aSystem(aCompartment->system);
  FOR_ALL(System::Variables, aSystem->getVariables())
    {
      Variable* aVariable(i->second);
      if(aVariable->getID() == "LIPID" || aVariable->getID() == "VACANT")
        {
          if(aVariable->getValue())
            {
              aCompartment->isEnclosed = true;
            }
          else
            {
              aCompartment->isEnclosed = false;
            }
          //Set the number of lipid/vacant molecules to be always 0 because
          //when we populate lattice we shouldn't create more lipid/vacant
          //molecules than the ones already created for the compartment:
          aVariable->setValue(0);
          Species* aSpecies(addSpecies(aVariable));
          aCompartment->vacantID = aSpecies->getID();
          if(aVariable->getID() == "LIPID")
            { 
              aSpecies->setIsLipid();
            }
          else
            {
              aSpecies->setIsVacant();
            }
        }
      std::vector<Species*>::iterator j(variable2species(aVariable));
      if(j != theSpecies.end())
        {
          aCompartment->species.push_back(*j);
          (*j)->setCompartment(aCompartment);
          if(!aCompartment->isSurface)
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
  Compartment* aRootCompartment(theCompartments[0]);
  theHCPk = theNormalizedVoxelRadius/sqrt(3); 
  theHCPh = theNormalizedVoxelRadius*sqrt(8.0/3);
  theHCPl = theNormalizedVoxelRadius*sqrt(3);
  if(aRootCompartment->shape == SPHERICAL || aRootCompartment->shape == ROD ||
     aRootCompartment->shape == ELLIPSOID)
    {
      switch(LatticeType)
        {
        case HCP_LATTICE: 
          theCenterPoint.z = aRootCompartment->lengthZ/2+4*
            theNormalizedVoxelRadius; //row
          theCenterPoint.y = aRootCompartment->lengthY/2+2*theHCPl; //layer
          theCenterPoint.x = aRootCompartment->lengthX/2+2*theHCPh; //column
          break;
        case CUBIC_LATTICE:
          theCenterPoint.z = aRootCompartment->lengthZ/2+4*
            theNormalizedVoxelRadius; //row
          theCenterPoint.y = aRootCompartment->lengthY/2+4*
            theNormalizedVoxelRadius; //layer
          theCenterPoint.x = aRootCompartment->lengthX/2+4*
            theNormalizedVoxelRadius; //column
          break;
        }
    }
  else if(aRootCompartment->shape == CUBIC || aRootCompartment->shape == CUBOID)
    {
      //We do not give any leeway space between the simulation boundary
      //and the cell boundary if it is CUBIC or CUBOID to support
      //periodic boundary conditions:
      theCenterPoint.z = aRootCompartment->lengthZ/2; //row
      theCenterPoint.y = aRootCompartment->lengthY/2; //layer
      theCenterPoint.x = aRootCompartment->lengthX/2; //column
    }
  aRootCompartment->centerPoint = theCenterPoint; 
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
  //For the CUBIC and CUBOID cell shapes, we need to readjust the size of
  //row, layer and column according to the boundary condition of its surfaces
  //to reflect the correct volume. This is because periodic boundary will
  //[consume a layer of the surface voxels:
  if(aRootCompartment->shape == CUBIC || aRootCompartment->shape == CUBOID)
    {
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
  for(unsigned int i(0); i != theCompartments.size(); ++i)
    {
      Compartment* aCompartment(theCompartments[i]); 
      if(aCompartment->isSurface)
        {
          aCompartment->actualArea =  (72*pow(VoxelRadius,2))*
            aCompartment->coords.size()/(6*pow(2,0.5)+4*pow(3,0.5)+
                                         3*pow(6, 0.5));
        }
      else
        { 
          int voxelCnt(aCompartment->coords.size());
          for(unsigned int j(0); j != aCompartment->allSubs.size(); ++j)
            {
              voxelCnt += aCompartment->allSubs[j]->coords.size();
            }
          aCompartment->actualVolume = (4*pow(2,0.5)*pow(VoxelRadius,3))*
            voxelCnt;
        }
    }
}

void SpatiocyteStepper::printSimulationParameters()
{
  std::cout << std::endl;
  for(unsigned int i(0); i != theSpecies.size()-1; ++i)
    {
      std::cout << "id:" << i << " " << 
        theSpecies[i]->getVariable()->getFullID().asString() << std::endl;
    }
  std::cout << "id:" << theSpecies.size()-1  << " NULL" << std::endl; 
  std::cout << "Voxel radius, r_v:" << VoxelRadius << " m" << std::endl;
  std::cout << "Simulation height:" << theCenterPoint.y*2*VoxelRadius*2 <<
    " m" << std::endl;
  std::cout << "Simulation width:" << theCenterPoint.z*2*VoxelRadius*2 << 
    " m" << std::endl;
  std::cout << "Simulation length:" << theCenterPoint.x*2*VoxelRadius*2 <<
    " m" << std::endl;
  std::cout << "Row size:" << theRowSize << std::endl;
  std::cout << "Layer size:" << theLayerSize << std::endl;
  std::cout << "Column size:" << theColSize << std::endl;
  std::cout << "Total allocated voxels:" << 
    theRowSize*theLayerSize*theColSize << std::endl;
  for(unsigned int i(0); i != theCompartments.size(); ++i)
    {
      Compartment* aCompartment(theCompartments[i]);
      double aSpecVolume(aCompartment->specVolume);
      double aSpecArea(aCompartment->specArea);
      double anActualVolume(aCompartment->actualVolume);
      double anActualArea(aCompartment->actualArea);
      switch(aCompartment->shape)
        {
        case SPHERICAL:
          std::cout << "Spherical (radius=" << 
            pow(3*aSpecVolume/(4*M_PI), 1.0/3) << "m) ";
          break;
        case ROD:
          std::cout << "Rod (radius=" << aCompartment->lengthY*VoxelRadius << 
            "m, cylinder length=" <<
            (aCompartment->eastPoint.x-aCompartment->westPoint.y)*
            VoxelRadius*2 << "m) ";
          break;
        case CUBIC:
          std::cout << "Cubic ";
          break;
        case CUBOID:
          std::cout << "Cuboid ";
          break;
        case ELLIPSOID:
          std::cout << "Ellipsoid ";
          break;
        }
      std::cout << aCompartment->system->getFullID().asString();
      if(aCompartment->isSurface)
        {
          std::cout << " Surface Compartment:" << std::endl;
          std::cout << "  [" << int(aSpecArea*(6*sqrt(2)+4*sqrt(3)+3*sqrt(6))/
                              (72*VoxelRadius*VoxelRadius)) << 
            "] Specified surface voxels {n_s = S_specified*"
            << "(6*2^0.5+4*3^0.5+3*6^0.5)/(72*r_v^2}" << std::endl;
          std::cout << "  [" << aCompartment->coords.size() <<
            "] Actual surface voxels {n_s}" << std::endl;
          std::cout << "  [" << aSpecArea << " m^2] Specified surface area " <<
            "{S_specified}" << std::endl;
          std::cout << "  [" << anActualArea << " m^2] Actual surface area " <<
            "{S = (72*r_v^2)*n_s/(6*2^0.5+4*3^0.5+3*6^0.5)}" << std::endl;
        }
      else
        {
          std::cout << " Volume Compartment:" << std::endl;
          int voxelCnt(aCompartment->coords.size());
          for(unsigned int j(0); j != aCompartment->allSubs.size(); ++j)
            {
              voxelCnt += aCompartment->allSubs[j]->coords.size();
            }
          std::cout << "  [" << int(aSpecVolume/(4*sqrt(2)*pow(VoxelRadius, 3))) << 
            "] Specified volume voxels {n_v = V_specified/(4*2^0.5*r_v^3)}" <<
          std::endl;  
          std::cout << "  [" << voxelCnt << "] Actual volume voxels {n_v}"  << std::endl;
          std::cout << "  [" << aSpecVolume << " m^3] Specified volume {V_specified}"
            << std::endl; 
          std::cout << "  [" << anActualVolume << " m^3] Actual volume " <<
            "{V = (4*2^0.5*r_v^3)*n_v}" << std::endl; 
        }
    }
  std::cout << std::endl;
  printProcessParameters();
}

void SpatiocyteStepper::readjustSurfaceBoundarySizes()
{
  Compartment* aRootCompartment(theCompartments[0]);
  //[XY, XZ, YZ]PLANE: the boundary type of the surface when 
  //the shape of the root compartment is CUBIC or CUBOID.
  //Boundary type can be either PERIODIC or REFLECTIVE.
  //Increase the size of [row,layer,col] by one voxel and make them odd sized
  //if the system uses periodic boundary conditions.
  if(aRootCompartment->yzPlane == PERIODIC)
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
  if(aRootCompartment->xzPlane == PERIODIC)
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
  if(aRootCompartment->xyPlane == PERIODIC)
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

void SpatiocyteStepper::constructLattice()
{ 
  Compartment* aRootCompartment(theCompartments[0]);
  unsigned int aSize(theRowSize*theLayerSize*theColSize);
  unsigned int a(0);
  unsigned int b(theStartCoord);
  unsigned short rootID(aRootCompartment->vacantID);
  for(std::vector<Voxel>::iterator i(theLattice.begin()); a != aSize; ++i, ++a, ++b)
    { 
      unsigned int aCol(a/(theRowSize*theLayerSize)); 
      unsigned int aLayer((a%(theRowSize*theLayerSize))/theRowSize); 
      unsigned int aRow((a%(theRowSize*theLayerSize))%theRowSize); 
      (*i).coord = b; 
      if(aRootCompartment->shape == CUBIC || 
         aRootCompartment->shape == CUBOID ||
         isInsideCoord(b, aRootCompartment, 0))
        {
          //By default, the voxel is vacant and we set it to the root id:
          (*i).id = rootID;
          for(unsigned int j(0); j != ADJOINING_VOXEL_SIZE; ++j)
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
          //Concatenate the some of the null voxels close to the surface:
          if(isInsideCoord(b, aRootCompartment, -4))
            {
              concatenateVoxel(&(*i), aRow, aLayer, aCol);
            }
        }
    }
  if(aRootCompartment->shape == CUBIC || aRootCompartment->shape == CUBOID)
    {
      concatenatePeriodicSurfaces();
    }
}

void SpatiocyteStepper::setPeriodicEdge()
{
  isPeriodicEdge = true;
}

bool SpatiocyteStepper::isPeriodicEdgeCoord(unsigned int aCoord,
                                            Compartment* aCompartment)
{
  unsigned int aRow;
  unsigned int aLayer;
  unsigned int aCol;
  coord2global(aCoord, &aRow, &aLayer, &aCol);
  if(aCompartment->system->getSuperSystem()->isRootSystem() &&
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

bool SpatiocyteStepper::isRemovableEdgeCoord(unsigned int aCoord,
                                             Compartment* aCompartment)
{
  unsigned int aRow;
  unsigned int aLayer;
  unsigned int aCol;
  coord2global(aCoord, &aRow, &aLayer, &aCol);
  if((aCompartment->xyPlane == REMOVE_UPPER && aRow == theRowSize-2) ||
     (aCompartment->xyPlane == REMOVE_LOWER && aRow == 1) ||
     (aCompartment->xyPlane == REMOVE_BOTH &&
      (aRow == 1 || aRow == theRowSize-2)) ||
     (aCompartment->xzPlane == REMOVE_UPPER && aLayer == theLayerSize-2) ||
     (aCompartment->xzPlane == REMOVE_LOWER && aLayer == 1) ||
     (aCompartment->xzPlane == REMOVE_BOTH &&
      (aLayer == 1 || aLayer == theLayerSize-2)) ||
     (aCompartment->yzPlane == REMOVE_UPPER && aCol == theColSize-2) ||
     (aCompartment->yzPlane == REMOVE_LOWER && aCol == 1) ||
     (aCompartment->yzPlane == REMOVE_BOTH &&
      (aCol == 1 || aCol == theColSize-2)))
    {
      return true;
    }
  return false;
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

void SpatiocyteStepper::setCompartmentProperties(Compartment* aCompartment)
{
  System* aSystem(aCompartment->system);
  double aRadius(0);
  switch(aCompartment->shape)
    {
    case SPHERICAL:
      if(!aCompartment->specVolume)
        {
          THROW_EXCEPTION(NotFound, "Property SIZE of the Sphere compartment "
                          + aSystem->getFullID().asString() + " not defined" );
        }
      aRadius = pow(3*aCompartment->specVolume/(4*M_PI), 1.0/3);
      aCompartment->lengthX = 2*aRadius;
      aCompartment->lengthY = aCompartment->lengthX;
      aCompartment->lengthZ = aCompartment->lengthX;
      aCompartment->specArea = 4*M_PI*aRadius*aRadius;
      break;
    case ROD:
      if(!aCompartment->specVolume)
        {
          THROW_EXCEPTION(NotFound, "Property SIZE of the Rod compartment "
                          + aSystem->getFullID().asString() + " not defined." );
        }
      if(!aCompartment->lengthY)
        {
          THROW_EXCEPTION(NotFound, "Property LENGTHY of the Rod compartment "
                          + aSystem->getFullID().asString() + ", which "
                          + "specifies its diameter, not defined." );
        }
      aRadius = aCompartment->lengthY/2;
      aCompartment->lengthX = aCompartment->specVolume/
        (M_PI*aRadius*aRadius)-(4*aRadius/3)+(2*aRadius);
      aCompartment->lengthZ = aCompartment->lengthY;
      aCompartment->specArea = 4*M_PI*aRadius*aRadius+
        2*M_PI*aRadius*(aCompartment->lengthX-2*aRadius);
      break;
    case CUBIC:
      if(!aCompartment->specVolume)
        {
          THROW_EXCEPTION(NotFound, "Property SIZE of the Cubic compartment "
                          + aSystem->getFullID().asString() + " not defined.");
        }
      aCompartment->lengthX = pow(aCompartment->specVolume, 1.0/3);
      aCompartment->lengthY = aCompartment->lengthX;
      aCompartment->lengthZ = aCompartment->lengthX;
      aCompartment->specArea = getCuboidSpecArea(aCompartment);
      break;
    case CUBOID:
      if(!aCompartment->lengthX)
        {
          THROW_EXCEPTION(NotFound,
                          "Property LENGTHX of the Cuboid compartment " +
                          aSystem->getFullID().asString() + " is not defined.");
        }
      if(!aCompartment->lengthY)
        {
          THROW_EXCEPTION(NotFound,
                          "Property LENGTHY of the Cuboid compartment " +
                          aSystem->getFullID().asString() + " is not defined.");
        }
      if(!aCompartment->lengthZ)
        {
          THROW_EXCEPTION(NotFound,
                          "Property LENGTHZ of the Cuboid compartment " +
                          aSystem->getFullID().asString() + " is not defined.");
        }
      aCompartment->specVolume = aCompartment->lengthX*aCompartment->lengthY*
        aCompartment->lengthZ;
      aCompartment->specArea = getCuboidSpecArea(aCompartment);
      break;
    case ELLIPSOID:
      if(!aCompartment->lengthX)
        {
          THROW_EXCEPTION(NotFound,
                          "Property LENGTHX of the Ellipsoid compartment " +
                          aSystem->getFullID().asString() + " is not defined.");
        }
      if(!aCompartment->lengthY)
        {
          THROW_EXCEPTION(NotFound,
                          "Property LENGTHY of the Ellipsoid compartment " +
                          aSystem->getFullID().asString() + " is not defined.");
        }
      if(!aCompartment->lengthZ)
        {
          THROW_EXCEPTION(NotFound,
                          "Property LENGTHZ of the Ellipsoid compartment " +
                          aSystem->getFullID().asString() + " is not defined.");
        }
      aCompartment->specVolume = 4*M_PI*aCompartment->lengthX*
        aCompartment->lengthY* aCompartment->lengthZ/24;
      aCompartment->specArea = 4*M_PI*
        pow((pow(aCompartment->lengthX/2, 1.6075)*
             pow(aCompartment->lengthY/2, 1.6075)+
             pow(aCompartment->lengthX/2, 1.6075)*
             pow(aCompartment->lengthZ/2, 1.6075)+
             pow(aCompartment->lengthY/2, 1.6075)*
             pow(aCompartment->lengthZ/2, 1.6075))/ 3, 1/1.6075); 
      break;
    }
  aCompartment->lengthX /= VoxelRadius*2;
  aCompartment->lengthY /= VoxelRadius*2;
  aCompartment->lengthZ /= VoxelRadius*2;
}

void SpatiocyteStepper::setCompartmentCenterPoint(Compartment* aCompartment)
{
  System* aSystem(aCompartment->system);
  System* aSuperSystem(aSystem->getSuperSystem());
  if(aCompartment->isSurface)
    {
      aSystem = aCompartment->system->getSuperSystem();
      aSuperSystem = aSystem->getSuperSystem();
    }
  Compartment* aSuperCompartment(system2compartment(aSuperSystem));
  //The center with reference to the immediate super system:
  aCompartment->centerPoint = aSuperCompartment->centerPoint;
  aCompartment->centerPoint.x += 
    aCompartment->originX*aSuperCompartment->lengthX/2;
  aCompartment->centerPoint.y += 
    aCompartment->originY*aSuperCompartment->lengthY/2;
  aCompartment->centerPoint.z += 
    aCompartment->originZ*aSuperCompartment->lengthZ/2;
  //The center of the west hemisphere of the rod:
  aCompartment->westPoint.x = aCompartment->centerPoint.x-
    aCompartment->lengthX/2+aCompartment->lengthY/2;
  aCompartment->westPoint.y = aCompartment->centerPoint.y;
  aCompartment->westPoint.z = aCompartment->centerPoint.z;
  //The center of the east hemisphere of the rod
  aCompartment->eastPoint.x = aCompartment->centerPoint.x+
    aCompartment->lengthX/2-aCompartment->lengthY/2;
  aCompartment->eastPoint.y = aCompartment->centerPoint.y;
  aCompartment->eastPoint.z = aCompartment->centerPoint.z;
  if(aCompartment->shape == CUBOID)
    {
      //The center of the west hemisphere of the rod:
      aCompartment->westPoint.x = aCompartment->centerPoint.x-
        aCompartment->lengthX/2;
      aCompartment->westPoint.y = aCompartment->centerPoint.y-
        aCompartment->lengthY/2;
      aCompartment->westPoint.z = aCompartment->centerPoint.z-
        aCompartment->lengthZ/2;
      //The center of the east hemisphere of the rod
      aCompartment->eastPoint.x = aCompartment->centerPoint.x+
        aCompartment->lengthX/2;
      aCompartment->eastPoint.y = aCompartment->centerPoint.y+
        aCompartment->lengthY/2;
      aCompartment->eastPoint.z = aCompartment->centerPoint.z+
        aCompartment->lengthZ/2;
    }
}

double SpatiocyteStepper::getCuboidSpecArea(Compartment* aCompartment)
{
  double anArea(0);
  if(aCompartment->xyPlane == UNIPERIODIC || 
     aCompartment->xyPlane == REMOVE_UPPER ||
     aCompartment->xyPlane == REMOVE_LOWER)
    { 
      anArea += aCompartment->lengthX*aCompartment->lengthY; 
    }
  else if(aCompartment->xyPlane == REFLECTIVE)
    {
      anArea += 2*aCompartment->lengthX*aCompartment->lengthY; 
    }
  if(aCompartment->xzPlane == UNIPERIODIC || 
     aCompartment->xzPlane == REMOVE_UPPER ||
     aCompartment->xzPlane == REMOVE_LOWER)
    { 
      anArea += aCompartment->lengthX*aCompartment->lengthZ; 
    }
  else if(aCompartment->xzPlane == REFLECTIVE)
    {
      anArea += 2*aCompartment->lengthX*aCompartment->lengthZ; 
    }
  if(aCompartment->yzPlane == UNIPERIODIC || 
     aCompartment->yzPlane == REMOVE_UPPER ||
     aCompartment->yzPlane == REMOVE_LOWER)
    { 
      anArea += aCompartment->lengthY*aCompartment->lengthZ; 
    }
  else if(aCompartment->yzPlane == REFLECTIVE)
    {
      anArea += 2*aCompartment->lengthY*aCompartment->lengthZ; 
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
  unsigned int aGlobalCol((unsigned int)(aPoint.x/theHCPh));
  unsigned int aGlobalLayer((unsigned int)((aPoint.y-(aGlobalCol%2)*theHCPk)/
                                           theHCPl));
  unsigned int aGlobalRow((unsigned int)((aPoint.z-
          ((aGlobalLayer+aGlobalCol)%2)*
          theNormalizedVoxelRadius)/(2*theNormalizedVoxelRadius)));
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
  //Specific for HCP lattice:{
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
  //}
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
  //Specific for HCP lattice:{
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
  //}
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
  Compartment* aRootCompartment(theCompartments[0]);
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
      if(aRootCompartment->xyPlane == UNIPERIODIC)
        { 
          replaceUniVoxel(aSrcVoxel, aDestVoxel);
        }
      else if(aRootCompartment->xyPlane == PERIODIC)
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

      if(aRootCompartment->xzPlane == UNIPERIODIC)
        { 
          replaceUniVoxel(aSrcVoxel, aDestVoxel);
        }
      else if(aRootCompartment->xzPlane == PERIODIC)
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
      if(aRootCompartment->yzPlane == UNIPERIODIC)
        { 
          replaceUniVoxel(aSrcVoxel, aDestVoxel);
        }
      else if(aRootCompartment->yzPlane == PERIODIC)
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
      for(unsigned int j(0); j!=ADJOINING_VOXEL_SIZE; ++j)
        {
          if(aSrcVoxel->adjoiningVoxels[j] == aSrcVoxel &&
              aDestVoxel->adjoiningVoxels[j] != aDestVoxel)
            {
              aSrcVoxel->adjoiningVoxels[j] = aDestVoxel->adjoiningVoxels[j];
              for(unsigned int k(0); k!=ADJOINING_VOXEL_SIZE; ++k)
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
      for(unsigned int j(0); j!=ADJOINING_VOXEL_SIZE; ++j)
        {
          if(aSrcVoxel->adjoiningVoxels[j] == aSrcVoxel)
            {
              for(unsigned int k(0); k!=ADJOINING_VOXEL_SIZE; ++k)
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
          gsl_ran_shuffle(getRng(), (*i).adjoiningVoxels, ADJOINING_VOXEL_SIZE,
                          sizeof(Voxel*));
        }
    }
}

void SpatiocyteStepper::setSurfaceVoxelProperties()
{
  for(std::vector<Compartment*>::iterator i(theCompartments.begin());
      i != theCompartments.end(); ++i)
    {
      if((*i)->isSurface)
        {
          setReactiveCompartments(*i);
          removePeriodicEdgeVoxels(*i);
          removeSurfaces(*i);
          for(std::vector<unsigned int>::iterator j((*i)->coords.begin());
              j != (*i)->coords.end(); ++j)
            {
              Voxel* aVoxel(&theLattice[*j]);
              optimizeSurfaceVoxel(aVoxel, *i);
              setSurfaceSubunit(aVoxel, *i);
            }
        }
    }
}

void SpatiocyteStepper::removePeriodicEdgeVoxels(Compartment* aCompartment)
{ 
  if(isPeriodicEdge)
    {
      std::vector<unsigned int> coords;
      for(std::vector<unsigned int>::iterator j(aCompartment->coords.begin());
          j != aCompartment->coords.end(); ++j)
        {
          Voxel* aVoxel(&theLattice[*j]);
          if(isPeriodicEdgeCoord(aVoxel->coord, aCompartment))
            {
              aVoxel->id = theNullID;
            }
          else
            {
              coords.push_back(*j);
            }
        }
      aCompartment->coords = coords;
    }
}

void SpatiocyteStepper::removeSurfaces(Compartment* aCompartment)
{ 
  std::vector<unsigned int> coords;
  for(std::vector<unsigned int>::iterator j(aCompartment->coords.begin());
      j != aCompartment->coords.end(); ++j)
    {
      Voxel* aVoxel(&theLattice[*j]);
      if(isRemovableEdgeCoord(aVoxel->coord, aCompartment))
        {
          aVoxel->id = theNullID;
        }
      else
        {
          coords.push_back(*j);
        }
    }
  aCompartment->coords = coords;
}

void SpatiocyteStepper::optimizeSurfaceVoxel(Voxel* aVoxel,
                                             Compartment* aCompartment)
{
  unsigned short surfaceID(aCompartment->vacantID);
  aVoxel->surfaceVoxels = new std::vector<std::vector<Voxel*> >;
  aVoxel->surfaceVoxels->resize(4);
  std::vector<Voxel*>& immediateSurface((*aVoxel->surfaceVoxels)[IMMEDIATE]);
  std::vector<Voxel*>& extendedSurface((*aVoxel->surfaceVoxels)[EXTENDED]);
  std::vector<Voxel*>& innerVolume((*aVoxel->surfaceVoxels)[INNER]);
  std::vector<Voxel*>& outerVolume((*aVoxel->surfaceVoxels)[OUTER]);
  std::vector<std::vector<Voxel*> > sharedVoxelsList;
  Voxel** forward(aVoxel->adjoiningVoxels);
  Voxel** reverse(forward+ADJOINING_VOXEL_SIZE);
  std::vector<Voxel*> adjoiningCopy;
  for(int k(0); k != ADJOINING_VOXEL_SIZE; ++k)
    {
      adjoiningCopy.push_back(forward[k]);
    }
  //Separate adjoining surface voxels and adjoining volume voxels.
  //Put the adjoining surface voxels at the beginning of the
  //adjoiningVoxels list while the volume voxels are put at the std::endl:
  for(std::vector<Voxel*>::iterator l(adjoiningCopy.begin());
      l != adjoiningCopy.end(); ++l)
    {
      if((*l) != aVoxel && ((*l)->id == surfaceID || 
                isReactiveCompartment(aCompartment, id2compartment((*l)->id))))
        {
          (*forward) = (*l);
          ++forward;
          //immediateSurface contains all adjoining surface voxels except the 
          //source voxel, aVoxel:
          immediateSurface.push_back(*l);
          for(int m(0); m != ADJOINING_VOXEL_SIZE; ++m)
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
          //surface (i.e., the parent volume) compartment, it must be an
          //outer volume voxel:
          if(!isInsideCoord((*l)->coord, aCompartment, 0))
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

void SpatiocyteStepper::setReactiveCompartments(Compartment* aCompartment)
{
  FOR_ALL(System::Variables, aCompartment->system->getVariables())
    {
      Variable* aVariable(i->second);
      if(aVariable->getID() == "REACTIVE")
        {
          String aStringID(aVariable->getName()); 
          aStringID = "System:" + aStringID;
          FullID aFullID(aStringID);
          System* aSystem(static_cast<System*>(getModel()->getEntity(aFullID)));
          aCompartment->reactiveComps.push_back(system2compartment(aSystem));
        }
    }
}

bool SpatiocyteStepper::isReactiveCompartment(Compartment* sourceCompartment,
                                              Compartment* targetCompartment)
{
  for(unsigned int i(0); i != sourceCompartment->reactiveComps.size(); ++i) 
    {
      if(sourceCompartment->reactiveComps[i] == targetCompartment)
        {
          return true;
        }
    }
  return false;
}

Species* SpatiocyteStepper::id2species(unsigned short id)
{
  return theSpecies[id];
}

Compartment* SpatiocyteStepper::id2compartment(unsigned short id)
{
  return theSpecies[id]->getCompartment();
}

Voxel* SpatiocyteStepper::coord2voxel(unsigned int aCoord)
{
  return &theLattice[aCoord];
}

Compartment* SpatiocyteStepper::system2compartment(System* aSystem)
{
  for(unsigned int i(0); i != theCompartments.size(); ++i)
    {
      if(theCompartments[i]->system == aSystem)
        {
          return theCompartments[i];
        }
    }
  return NULL;
}

void SpatiocyteStepper::setSurfaceSubunit(Voxel* aVoxel,
                                          Compartment* aCompartment)
{
  aVoxel->subunit = new Subunit;
  aVoxel->subunit->voxel = aVoxel;
  Point& aPoint(aVoxel->subunit->surfacePoint);
  aPoint = coord2point(aVoxel->coord);
  double aRadius(aCompartment->lengthY/2);
  Point aCenterPoint(aCompartment->centerPoint);
  if(aPoint.x < aCompartment->westPoint.x)
    {
      aCenterPoint.x = aCompartment->westPoint.x;
    }
  else if(aPoint.x > aCompartment->eastPoint.x )
    {
      aCenterPoint.x = aCompartment->eastPoint.x;
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

void SpatiocyteStepper::compartmentalizeLattice() 
{
  for(std::vector<Voxel>::iterator i(theLattice.begin()); i != theLattice.end(); ++i)
    {
      if((*i).id != theNullID)
        { 
          compartmentalizeVoxel(&(*i), theCompartments[0]);
        }
    }
}


bool SpatiocyteStepper::compartmentalizeVoxel(Voxel* aVoxel,
                                              Compartment* aCompartment)
{
  if(!aCompartment->isSurface)
    {
      if(aCompartment->system->isRootSystem() ||
         isInsideCoord(aVoxel->coord, aCompartment, 0))
        {
          for(unsigned int i(0); i != aCompartment->immediateSubs.size(); ++i)
            {
              if(compartmentalizeVoxel(aVoxel, aCompartment->immediateSubs[i]))
                {
                  return true;
                }
            }
          if(aCompartment->surfaceSub && 
             (aCompartment->system->isRootSystem() ||
              aCompartment->surfaceSub->isEnclosed))
            {
              if(isEnclosedSurfaceVoxel(aVoxel, aCompartment))
                {
                  aVoxel->id = aCompartment->surfaceSub->vacantID;
                  aCompartment->surfaceSub->coords.push_back(
                                                 aVoxel->coord-theStartCoord);
                  return true;
                }
            }
          aVoxel->id = aCompartment->vacantID;
          aCompartment->coords.push_back(aVoxel->coord-theStartCoord);
          return true;
        }
      if(aCompartment->surfaceSub)
        {
          if(isSurfaceVoxel(aVoxel, aCompartment))
            {
              aVoxel->id = aCompartment->surfaceSub->vacantID;
              aCompartment->surfaceSub->coords.push_back(
                                         aVoxel->coord-theStartCoord);
              return true;
            }
        }
    }
  return false;
}

bool SpatiocyteStepper::isSurfaceVoxel(Voxel* aVoxel, Compartment* aCompartment)
{
  for(unsigned int i(0); i != ADJOINING_VOXEL_SIZE; ++i)
    {
      if(isInsideCoord(aVoxel->adjoiningVoxels[i]->coord, aCompartment, 0))
        {
          return true;
        }
    }
  return false;
}


bool SpatiocyteStepper::isEnclosedSurfaceVoxel(Voxel* aVoxel,
                                               Compartment* aCompartment)
{
  Compartment* aSuperCompartment(system2compartment(
                          aCompartment->system->getSuperSystem()));
  for(unsigned int i(0); i != ADJOINING_VOXEL_SIZE; ++i)
    {
      if(aVoxel->adjoiningVoxels[i]->id == theNullID ||
         aVoxel->adjoiningVoxels[i] == aVoxel ||
         !isInsideCoord(aVoxel->adjoiningVoxels[i]->coord,
                        aSuperCompartment, 0))
        {
          return true;
        }
    }
  return false;
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
                                      Compartment* aCompartment, double delta)
{
  Point aPoint(coord2point(aCoord));
  Point aCenterPoint(aCompartment->centerPoint);
  aPoint.x = aPoint.x - aCenterPoint.x;
  aPoint.y = aPoint.y - aCenterPoint.y;
  aPoint.z = aPoint.z - aCenterPoint.z;
  rotateX(aCompartment->rotateX, &aPoint);
  rotateY(aCompartment->rotateY, &aPoint);
  rotateZ(aCompartment->rotateZ, &aPoint);
  aPoint.x = aPoint.x + aCenterPoint.x;
  aPoint.y = aPoint.y + aCenterPoint.y;
  aPoint.z = aPoint.z + aCenterPoint.z;
  double aRadius(aCompartment->lengthY/2+theNormalizedVoxelRadius-delta);
  switch(aCompartment->shape)
    {
    case ELLIPSOID:
    case SPHERICAL:
      //If the distance between the voxel and the center point is less than 
      //or equal to radius-2, then the voxel cannot be a surface voxel:
      if(pow(aPoint.x-aCenterPoint.x, 2)/
         pow((aCompartment->lengthX-delta)/2, 2)+ 
         pow(aPoint.y-aCenterPoint.y, 2)/
         pow((aCompartment->lengthY-delta)/2, 2)+ 
         pow(aPoint.z-aCenterPoint.z, 2)/
         pow((aCompartment->lengthZ-delta)/2, 2) <= 1)
        {
          return true;
        }
      break;
    case ROD: 
      //The axial point of the cylindrical portion of the rod:
      aCenterPoint.x = aPoint.x;
      //If the distance between the voxel and the center point is less than 
      //or equal to the radius, then the voxel must be inside the compartment:
      if((aPoint.x >= aCompartment->westPoint.x &&
          aPoint.x <= aCompartment->eastPoint.x &&
          getDistance(&aPoint, &aCenterPoint) <= aRadius) ||
         (aPoint.x < aCompartment->westPoint.x &&
          getDistance(&aPoint, &aCompartment->westPoint) <= aRadius) ||
         (aPoint.x > aCompartment->eastPoint.x &&
          getDistance(&aPoint, &aCompartment->eastPoint) <= aRadius))
        { 
          return true;
        }
      break;
    case CUBIC:
      if(sqrt(pow(aPoint.x-aCenterPoint.x, 2)) <= aRadius && 
         sqrt(pow(aPoint.y-aCenterPoint.y, 2)) <= aRadius && 
         sqrt(pow(aPoint.z-aCenterPoint.z, 2)) <= aRadius)
        {
          return true;
        }
      break;
    case CUBOID:
      if(sqrt(pow(aPoint.x-aCenterPoint.x, 2)) <= 
         aCompartment->lengthX/2+theNormalizedVoxelRadius-delta &&
         sqrt(pow(aPoint.y-aCenterPoint.y, 2)) <= 
         aCompartment->lengthY/2+theNormalizedVoxelRadius-delta &&
         sqrt(pow(aPoint.z-aCenterPoint.z, 2)) <= 
         aCompartment->lengthZ/2+theNormalizedVoxelRadius-delta)
        {
          return true;
        }
      break;
    }
  return false;
}

void SpatiocyteStepper::populateCompartment(Compartment* aCompartment)
{
  unsigned int populationSize(0);
  unsigned int gaussianPopulationSize(0);
  //First, populate gaussian distributed molecules, if any:
  for(std::vector<Species*>::const_iterator i(aCompartment->species.begin());
      i != aCompartment->species.end(); ++i)
    {
      populationSize += (unsigned int)(*i)->getPopulateMoleculeSize();
      if((*i)->getIsGaussianPopulation())
        {
          gaussianPopulationSize +=
            (unsigned int)(*i)->getPopulateMoleculeSize();
          (*i)->populateCompartmentGaussian();
        }
    }
  //Second, populate remaining vacant voxels uniformly:
  //If there are many molecules to be populated we need to
  //systematically choose the vacant voxels randomly from a list
  //of available vacant voxels of the compartment:
  if(populationSize > gaussianPopulationSize)
    {
      unsigned int count(0);
      if(double(populationSize)/aCompartment->coords.size() > 0.2)
        {
          unsigned int* populateVoxels(new unsigned int[populationSize]);
          unsigned int availableVoxelSize(aCompartment->coords.size());
          unsigned int* availableVoxels(new unsigned int [availableVoxelSize]); 
          for(unsigned int i(0); i != availableVoxelSize; ++i)
            {
              availableVoxels[i] = i;
            }
          gsl_ran_choose(getRng(), populateVoxels, populationSize,
                     availableVoxels, availableVoxelSize, sizeof(unsigned int));
          //gsl_ran_choose arranges the position ascending, so we need
          //to shuffle the order of voxel positions:
          gsl_ran_shuffle(getRng(), populateVoxels, populationSize,
                          sizeof(unsigned int)); 
          for(std::vector<Species*>::const_iterator i(aCompartment->species.begin());
              i != aCompartment->species.end(); ++i)
            {
              if(!(*i)->getIsGaussianPopulation())
                {
                  (*i)->populateCompartmentUniform(populateVoxels, &count);
                }
            }
          delete[] populateVoxels;
          delete[] availableVoxels;
        }
      //Otherwise, we select a random voxel from the list of compartment
      //voxels and check if it is vacant before occupying it with a 
      //molecule, iteratively. This makes it much faster to populate
      //large compartments with small number of molecules.
      else
        {
          for(std::vector<Species*>::const_iterator i(aCompartment->species.begin());
              i != aCompartment->species.end(); ++i)
            {
              if(!(*i)->getIsGaussianPopulation())
                {
                  (*i)->populateCompartmentUniformSparse();
                }
            }
        }
    }
}


void SpatiocyteStepper::clearCompartment(Compartment* aCompartment)
{
  for(std::vector<Species*>::const_iterator i(aCompartment->species.begin());
      i != aCompartment->species.end(); ++i)
    {
      (*i)->removeMolecules();
    }
}


std::vector<Compartment*> const& SpatiocyteStepper::getCompartments() const
{
    return theCompartments;
}

