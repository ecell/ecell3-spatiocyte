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
#include <Model.hpp>
#include <System.hpp>
#include "SpatiocyteStepper.hpp"
#include "SpatiocyteSpecies.hpp"
#include "SpatiocyteProcess.hpp"
#include "ReactionProcess.hpp"

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
  cout << "2. setting up lattice properties..." << endl;
  setLatticeProperties(); 
  //We need a compartment tree to assign the voxels to each compartment
  //and get the available number of vacant voxels. The compartmentalized
  //vacant voxels are needed to randomly place molecules according to the
  //compartment:
  cout << "3. creating compartments..." << endl;
  registerCompartments();
  setCompartmentsProperties();
  //All species have been created at this point, we initialize them now:
  cout << "4. initializing species..." << endl;
  initSpecies();
  cout << "5. initializing processes the second time..." << endl;
  initProcessSecond();
  cout << "7. constructing lattice..." << endl;
  constructLattice();
  cout << "8. shuffling adjoining voxels..." << endl;
  shuffleAdjoiningVoxels();
  cout << "9. compartmentalizing lattice..." << endl;
  compartmentalizeLattice();
  cout << "10. setting up properties of surface voxels..." << endl;
  setSurfaceVoxelProperties();
  cout << "11. populating compartments with molecules..." << endl;
  populateCompartments();
  storeSimulationParameters();
  //checkSurfaceCompartment();
  cout << "12. initializing processes the third time..." << endl;
  initProcessThird();
  cout << "13. initializing the priority queue..." << endl;
  initPriorityQueue();
  printSimulationParameters();
  cout << "14. initializing processes the fourth time..." << endl;
  initProcessFourth();
  initProcessLastOnce();
  cout << "15. simulation is started..." << endl;
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
  vector<Species*>::iterator aSpeciesIter(variable2species(aVariable));
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
  vector<Species*>::iterator aSpeciesIter(variable2species(aVariable));
  if(aSpeciesIter == theSpecies.end())
    {
      return NULL;
    }
  return *aSpeciesIter;
}

vector<Species*> SpatiocyteStepper::getSpecies()
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
  //Specific for HCP lattice:{
  aPoint.y = (aCol%2)*theHCPk+theHCPl*aLayer;
  aPoint.z = aRow*2*theNormalizedVoxelRadius+
    ((aLayer+aCol)%2)*theNormalizedVoxelRadius;
  aPoint.x = aCol*theHCPh;
  //}
  return aPoint;
}


vector<Species*>::iterator
SpatiocyteStepper::variable2species(Variable* aVariable)
{
  for(vector<Species*>::iterator i(theSpecies.begin());
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
  vector<int> list;
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
      cout << "i:" << i << " ";
      if(theSpecies[i]->getVariable() != NULL)
        {
          cout << theSpecies[i]->getVariable()->getFullID().asString();
          if(theSpecies[i]->getIsVolume())
            {
              volumeCnt += list[i];
            }
          else
            {
              surfaceCnt += list[i];
            }
        }
      cout << " cnt:" << list[i] << endl;
    }
  cout << "total volume:" << volumeCnt << endl;
  cout << "total surface:" << surfaceCnt << endl;
  cout << "total volume+surface:" << surfaceCnt+volumeCnt << endl;
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
          cout << "size:" << (*i)->coords.size() << endl;
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
      cout << i << ": " << surfaceCnt[i] << endl;
    }
  cout << "total:" << total << endl;
  cout << "size:" << theSpecies.back()->size() << endl;
}
*/

void SpatiocyteStepper::initSpecies()
{
  for(vector<Species*>::iterator i(theSpecies.begin());
      i != theSpecies.end(); ++i)
    {
      (*i)->initialize(theSpecies.size());
    }
}

void SpatiocyteStepper::initProcessSecond()
{
  for(vector<Process*>::const_iterator i(theProcessVector.begin());
      i != theProcessVector.end(); ++i)
    {      
      SpatiocyteProcess* aProcess(reinterpret_cast<SpatiocyteProcess*>(*i));
      aProcess->initializeSecond();
    }
}

void SpatiocyteStepper::printProcessParameters()
{
  for(vector<Process*>::const_iterator i(theProcessVector.begin());
      i != theProcessVector.end(); ++i)
    {      
      SpatiocyteProcess* aProcess(reinterpret_cast<SpatiocyteProcess*>(*i));
      aProcess->printParameters();
    }
}

void SpatiocyteStepper::initProcessThird()
{
  for(vector<Process*>::const_iterator i(theProcessVector.begin());
      i != theProcessVector.end(); ++i)
    {      
      SpatiocyteProcess* aProcess(reinterpret_cast<SpatiocyteProcess*>(*i));
      aProcess->initializeThird();
    }
}

void SpatiocyteStepper::initProcessFourth()
{
  for(vector<Process*>::const_iterator i(theProcessVector.begin());
      i != theProcessVector.end(); ++i)
    {      
      SpatiocyteProcess* aProcess(reinterpret_cast<SpatiocyteProcess*>(*i));
      aProcess->initializeFourth();
    }
  setStepInterval(thePriorityQueue.getTop()->getTime()-getCurrentTime());
}

void SpatiocyteStepper::initProcessLastOnce()
{
  for(vector<Process*>::const_iterator i(theProcessVector.begin());
      i != theProcessVector.end(); ++i)
    {      
      SpatiocyteProcess* aProcess(reinterpret_cast<SpatiocyteProcess*>(*i));
      aProcess->initializeLastOnce();
    }
}

void SpatiocyteStepper::initPriorityQueue()
{
  const double aCurrentTime(getCurrentTime());
  thePriorityQueue.clear();
  for(vector<Process*>::const_iterator i(theProcessVector.begin());
      i != theProcessVector.end(); ++i)
    {      
      SpatiocyteProcess*
        aProcess(reinterpret_cast<SpatiocyteProcess*>(*i));
      aProcess->setTime(aCurrentTime+aProcess->getStepInterval());
      aProcess->setPriorityQueue(&thePriorityQueue);
      //The following processes are inserted in the PriorityQueue and
      //executed at simulation steps according to their execution times:
      String aClassName(aProcess->getPropertyInterface().getClassName());
      if(aClassName == "DiffusionProcess" ||
         aClassName == "IteratingLogProcess" ||
         aClassName == "MoleculePopulateProcess" ||
         aClassName == "CoordinateLogProcess" ||
         aClassName == "VisualizationLogProcess" ||
         aClassName == "FluorescentProteinImagingProcess" ||
         aClassName == "OscillationAnalysisProcess" ||
         aClassName == "PeriodicBoundaryDiffusionProcess" ||
         aClassName == "SpatiocyteNextReactionProcess" ||
         aClassName == "PolymerFragmentationProcess")
        {
          aProcess->setQueueID(thePriorityQueue.push(aProcess));
        }
      //The following processes never interrupt other Processes.
      //We exclude them here and set up the interrupt for the remaining
      //processes. All processes which interrupt other processes have
      //the ReactionProcess as the base class.
      if(aClassName != "DiffusionProcess" && 
         aClassName != "MoleculePopulateProcess" &&
         aClassName != "IteratingLogProcess" &&
         aClassName != "CoordinateLogProcess" &&
         aClassName != "VisualizationLogProcess" &&
         aClassName != "FluorescentProteinImagingProcess" &&
         aClassName != "OscillationAnalysisProcess" &&
         aClassName != "PeriodicBoundaryDiffusionProcess" && 
         aClassName != "PolymerizationParameterProcess")
        {
          ReactionProcess*
            aReactionProcess(reinterpret_cast<ReactionProcess*>(*i));
          aReactionProcess->setInterrupt(theProcessVector, *i);
        }
    } 
}

void SpatiocyteStepper::populateCompartments()
{
  for(vector<Compartment*>::const_iterator i(theCompartments.begin());
      i != theCompartments.end(); ++i)
    {
      populateCompartment(*i);
    }
}

void SpatiocyteStepper::clearCompartments()
{
  for(vector<Compartment*>::const_iterator i(theCompartments.begin());
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
  vector<Compartment*> allSubs;
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
  vector<Compartment*> compartments(theCompartments[0]->immediateSubs);
  while(!compartments.empty())
    {
      vector<Compartment*> subCompartments;
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

Compartment* SpatiocyteStepper::registerCompartment(System* aSystem,
                                            vector<Compartment*>* allSubs)
{ 
  //We execute this function to register the System, and its subsystems
  //recursively.
  Compartment* aCompartment(new Compartment);
  aCompartment->system = aSystem;
  aCompartment->surfaceSub = NULL;
  if(aSystem->getVariable("TYPE")->getValue() != SURFACE)
    {
      aCompartment->isSurface = false;
    }
  registerCompartmentSpecies(aCompartment);
  //Systems contains all the subsystems of a System.
  //For example /membrane is the subsystem of /:
  FOR_ALL(System::Systems, aSystem->getSystems())
    {
      Compartment* aSubCompartment(registerCompartment(i->second, allSubs)); 
      allSubs->push_back(aSubCompartment);
      aCompartment->immediateSubs.push_back(aSubCompartment);
      if(aSubCompartment->system->getVariable("TYPE")->getValue() ==
         SURFACE)
        {
          aSubCompartment->isSurface = true;
          aCompartment->surfaceSub = aSubCompartment;
        }
    }
  aCompartment->allSubs = *allSubs;
  return aCompartment;
}

void SpatiocyteStepper::setCompartmentsProperties()
{
  theCompartments[0]->centerPoint = theCenterPoint;
  for(unsigned int i(0); i != theCompartments.size(); ++i)
    {
      cout << theCompartments[i]->system->getFullID().asString() << endl;
      setCompartmentProperties(theCompartments[i]);
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
      vector<Species*>::iterator j(variable2species(aVariable));
      if(j != theSpecies.end())
        { 
          aCompartment->species.push_back(*j);
          (*j)->setCompartment(aCompartment);
          if(aSystem->getVariable("TYPE")->getValue() == VOLUME)
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
  //The root system is the cell:
  System* aRootSystem(getModel()->getRootSystem());
  //theCellShape: the shape of the cell, either SPHERICAL, ROD, CUBIC,
  //CUBOID or ELLIPSOID:
  theCellShape = (unsigned int)aRootSystem->getVariable("SHAPE")->getValue();
  double aCellLengthX(0);
  double aCellLengthY(0);
  double aCellLengthZ(0);
  double aVolume;
  double aCellRadius;
  switch(theCellShape)
    {
    case SPHERICAL:
      //Convert the value of SIZE in liters to meter^3 for the simulation
      //volume:
      aVolume = aRootSystem->getVariable("SIZE")->getValue()*1e-3;
      aCellLengthX = 2*pow(3*aVolume/(4*M_PI), 1.0/3);
      aCellLengthY = aCellLengthX;
      aCellLengthZ = aCellLengthX;
      break;
    case ROD:
      //Convert the value of SIZE in liters to meter^3 for the simulation
      //volume:
      aVolume = aRootSystem->getVariable("SIZE")->getValue()*1e-3;
      //Length of the cell from one end of the pole to the other
      //(i.e., including the radii of both hemispheres): 
      aCellLengthY = aRootSystem->getVariable("LENGTHY")->getValue();
      aCellRadius = aCellLengthY/2;
      aCellLengthX = aVolume/(M_PI*aCellRadius*aCellRadius)-
        (4*aCellRadius/3)+(2*aCellRadius);
      aCellLengthZ = aCellLengthY;
      break;
    case CUBIC:
      //Convert the value of SIZE in liters to meter^3 for the simulation
      //volume:
      aVolume = aRootSystem->getVariable("SIZE")->getValue()*1e-3;
      aCellLengthX = pow(aVolume, 1.0/3);
      aCellLengthY = aCellLengthX;
      aCellLengthZ = aCellLengthX;
      break;
    case CUBOID:
      aCellLengthX = aRootSystem->getVariable("LENGTHX")->getValue();
      aCellLengthY = aRootSystem->getVariable("LENGTHY")->getValue();
      aCellLengthZ = aRootSystem->getVariable("LENGTHZ")->getValue();
      break;
    case ELLIPSOID:
      aCellLengthX = aRootSystem->getVariable("LENGTHX")->getValue();
      aCellLengthY = aRootSystem->getVariable("LENGTHY")->getValue();
      aCellLengthZ = aRootSystem->getVariable("LENGTHZ")->getValue();
      break;
    }
  //VoxelRadius: the radius of a Spatiocyte voxel, r_v in the paper.
  //VoxelRadius is the actual radius of voxels which is set according to
  //the radius of simulated molecule. Here we divide theCellLengthXYZ
  //with VoxelRadius to normalize the voxel length
  //to 1 (i.e, theNormalizedVoxelRadius = 0.5) to avoid
  //dealing with very small values during calculations.
  aCellLengthX /= VoxelRadius*2;
  aCellLengthY /= VoxelRadius*2;
  aCellLengthZ /= VoxelRadius*2;
  //Specific properties of an hexagonal close-packed (HCP) lattice are
  //as follows:{
  theHCPk = theNormalizedVoxelRadius/sqrt(3); 
  theHCPh = theNormalizedVoxelRadius*sqrt(8.0/3);
  theHCPl = theNormalizedVoxelRadius*sqrt(3);
  //theCenterPoint: the center point of the entire simulation space:
  if(theCellShape == SPHERICAL || theCellShape == ROD ||
     theCellShape == ELLIPSOID)
    {
      theCenterPoint.z = aCellLengthZ/2+4*theNormalizedVoxelRadius; //row
      theCenterPoint.y = aCellLengthY/2+2*theHCPl; //layer
      theCenterPoint.x = aCellLengthX/2+2*theHCPh; //column
    }
  else if(theCellShape == CUBIC || theCellShape == CUBOID)
    {
      //We do not give any leeway space between the simulation boundary
      //and the cell boundary if it is CUBIC or CUBOID to support
      //periodic boundary conditions:
      theCenterPoint.z = aCellLengthZ/2; //row
      theCenterPoint.y = aCellLengthY/2; //layer
      theCenterPoint.x = aCellLengthX/2; //column
    }
  theRowSize = (unsigned int)rint((theCenterPoint.z)/
                                  (theNormalizedVoxelRadius));
  theLayerSize = (unsigned int)rint((theCenterPoint.y*2)/theHCPl);
  theColSize = (unsigned int)rint((theCenterPoint.x*2)/theHCPh);
  //}
  //For the CUBIC and CUBOID cell shapes, we need to readjust the size of
  //row, layer and column according to the boundary condition of its surfaces
  //to reflect the correct volume. This is because periodic boundary will
  //[consume a layer of the surface voxels:
  if(theCellShape == CUBIC || theCellShape == CUBOID)
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
  cout << endl;
  for(unsigned int i(0); i != theSpecies.size()-1; ++i)
    {
      cout << "id:" << i << " " << 
        theSpecies[i]->getVariable()->getFullID().asString() << endl;
    }
  cout << "id:" << theSpecies.size()-1  << " NULL" << endl; 
  cout << "Voxel radius, r_v:" << VoxelRadius << " m" << endl;
  cout << "Simulation height:" << theCenterPoint.y*2*VoxelRadius*2 <<
    " m" << endl;
  cout << "Simulation width:" << theCenterPoint.z*2*VoxelRadius*2 << 
    " m" << endl;
  cout << "Simulation length:" << theCenterPoint.x*2*VoxelRadius*2 <<
    " m" << endl;
  cout << "Row size:" << theRowSize << endl;
  cout << "Layer size:" << theLayerSize << endl;
  cout << "Column size:" << theColSize << endl;
  cout << "Total allocated voxels:" << 
    theRowSize*theLayerSize*theColSize << endl;
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
          cout << "Spherical (radius=" << 
            pow(3*aSpecVolume/(4*M_PI), 1.0/3) << "m) ";
          break;
        case ROD:
          cout << "Rod (radius=" << aCompartment->lengthY*VoxelRadius << 
            "m, cylinder length=" <<
            (aCompartment->eastPoint.x-aCompartment->westPoint.y)*
            VoxelRadius*2 << "m) ";
          break;
        case CUBIC:
          cout << "Cubic ";
          break;
        case CUBOID:
          cout << "Cuboid ";
          break;
        case ELLIPSOID:
          cout << "Ellipsoid ";
          break;
        }
      cout << aCompartment->system->getFullID().asString();
      if(aCompartment->isSurface)
        {
          cout << " Surface Compartment:" << endl;
          cout << "  [" << int(aSpecArea*(6*sqrt(2)+4*sqrt(3)+3*sqrt(6))/
                              (72*VoxelRadius*VoxelRadius)) << 
            "] Specified surface voxels {n_s = S_specified*"
            << "(6*2^0.5+4*3^0.5+3*6^0.5)/(72*r_v^2}" << endl;
          cout << "  [" << aCompartment->coords.size() <<
            "] Actual surface voxels {n_s}" << endl;
          cout << "  [" << aSpecArea << " m^2] Specified surface area " <<
            "{S_specified}" << endl;
          cout << "  [" << anActualArea << " m^2] Actual surface area " <<
            "{S = (72*r_v^2)*n_s/(6*2^0.5+4*3^0.5+3*6^0.5)}" << endl;
        }
      else
        {
          cout << " Volume Compartment:" << endl;
          int voxelCnt(aCompartment->coords.size());
          for(unsigned int j(0); j != aCompartment->allSubs.size(); ++j)
            {
              voxelCnt += aCompartment->allSubs[j]->coords.size();
            }
          cout << "  [" << int(aSpecVolume/(4*sqrt(2)*pow(VoxelRadius, 3))) << 
            "] Specified volume voxels {n_v = V_specified/(4*2^0.5*r_v^3)}" <<
          endl;  
          cout << "  [" << voxelCnt << "] Actual volume voxels {n_v}"  << endl;
          cout << "  [" << aSpecVolume << " m^3] Specified volume {V_specified}"
            << endl; 
          cout << "  [" << anActualVolume << " m^3] Actual volume " <<
            "{V = (4*2^0.5*r_v^3)*n_v}" << endl; 
        }
    }
  cout << endl;
  printProcessParameters();
}

void SpatiocyteStepper::readjustSurfaceBoundarySizes()
{
  //SURFACE[X,Y,Z]: the boundary type of the surface when theCellShape is CUBIC.
  //Boundary type can be either PERIODIC or REFLECTIVE.
  //Increase the size of [row,layer,col] by one voxel and make them odd sized
  //if the system uses periodic boundary conditions.
  if(getModel()->getRootSystem()->getVariable("SURFACEX")->getValue() ==
     PERIODIC)
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
  else
    {
      theColSize += 2;
    }
  if(getModel()->getRootSystem()->getVariable("SURFACEY")->getValue() ==
     PERIODIC)
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
  else
    {
      theLayerSize += 2;
    }
  if(getModel()->getRootSystem()->getVariable("SURFACEZ")->getValue() ==
     PERIODIC)
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
  else
    {
      theRowSize += 2;
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
  for(vector<Voxel>::iterator i(theLattice.begin()); a != aSize; ++i, ++a, ++b)
    { 
      unsigned int aCol(a/(theRowSize*theLayerSize)); 
      unsigned int aLayer((a%(theRowSize*theLayerSize))/theRowSize); 
      unsigned int aRow((a%(theRowSize*theLayerSize))%theRowSize); 
      (*i).coord = b; 
      if(theCellShape == CUBIC || theCellShape == CUBOID ||
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
  if(theCellShape == CUBIC || theCellShape == CUBOID)
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

void SpatiocyteStepper::setCompartmentProperties(Compartment* aCompartment)
{
  System* aSystem(aCompartment->system);
  System* aSuperSystem(aSystem->getSuperSystem());
  if(aCompartment->isSurface)
    {
      aSystem = aCompartment->system->getSuperSystem();
      aSuperSystem = aSystem->getSuperSystem();
    }
  Compartment* aSuperCompartment(system2compartment(aSuperSystem));
  aCompartment->shape = (unsigned int)aSystem->getVariable("SHAPE")->getValue();
  double anOriX(aSystem->getVariable("ORIGINX")->getValue());
  double anOriY(aSystem->getVariable("ORIGINY")->getValue());
  double anOriZ(aSystem->getVariable("ORIGINZ")->getValue());
  double aRadius;
  double anArea(0);
  switch(aCompartment->shape)
    {
    case SPHERICAL:
      aCompartment->specVolume = aSystem->getVariable("SIZE")->getValue()*1e-3;
      aRadius = pow(3*aCompartment->specVolume/(4*M_PI), 1.0/3);
      aCompartment->lengthX = 2*aRadius;
      aCompartment->lengthY = aCompartment->lengthX;
      aCompartment->lengthZ = aCompartment->lengthX;
      aCompartment->specArea = 4*M_PI*aRadius*aRadius;
      break;
    case ROD:
      aCompartment->specVolume = aSystem->getVariable("SIZE")->getValue()*1e-3;
      aCompartment->lengthY = aSystem->getVariable("LENGTHY")->getValue();
      aRadius = aCompartment->lengthY/2;
      aCompartment->lengthX = aCompartment->specVolume/
        (M_PI*aRadius*aRadius)-(4*aRadius/3)+(2*aRadius);
      aCompartment->lengthZ = aCompartment->lengthY;
      aCompartment->specArea = 4*M_PI*aRadius*aRadius+
        2*M_PI*aRadius*(aCompartment->lengthX-2*aRadius);
      break;
    case CUBIC:
      aCompartment->specVolume = aSystem->getVariable("SIZE")->getValue()*1e-3;
      aCompartment->lengthX = pow(aCompartment->specVolume, 1.0/3);
      aCompartment->lengthY = aCompartment->lengthX;
      aCompartment->lengthZ = aCompartment->lengthX;
      aCompartment->specArea = 6*aCompartment->lengthX*aCompartment->lengthX;
      break;
    case CUBOID:
      aCompartment->lengthX = aSystem->getVariable("LENGTHX")->getValue();
      aCompartment->lengthY = aSystem->getVariable("LENGTHY")->getValue();
      aCompartment->lengthZ = aSystem->getVariable("LENGTHZ")->getValue();
      aCompartment->specVolume = aCompartment->lengthX*aCompartment->lengthY*
        aCompartment->lengthZ;
      if(aSystem->getVariable("SURFACEZ")->getValue() ==
         UNIPERIODIC)
        { 
          anArea += aCompartment->lengthY*aCompartment->lengthX; 
        }
      else if(aSystem->getVariable("SURFACEZ")->getValue()
              != PERIODIC)
        { 
          anArea += 2*aCompartment->lengthY*aCompartment->lengthX; 
        }
      if(aSystem->getVariable("SURFACEX")->getValue() ==
         UNIPERIODIC)
        { 
          anArea += aCompartment->lengthY*aCompartment->lengthZ; 
        }
      else if(aSystem->getVariable("SURFACEX")->getValue()
              != PERIODIC)
        { 
          anArea += 2*aCompartment->lengthY*aCompartment->lengthZ; 
        }
      if(aSystem->getVariable("SURFACEY")->getValue() ==
         UNIPERIODIC)
        { 
          anArea += aCompartment->lengthX*aCompartment->lengthZ; 
        }
      else if(aSystem->getVariable("SURFACEY")->getValue()
              != PERIODIC)
        { 
          anArea += 2*aCompartment->lengthX*aCompartment->lengthZ; 
        }
      aCompartment->specArea = anArea;
      break;
    case ELLIPSOID:
      aCompartment->lengthX = aSystem->getVariable("LENGTHX")->getValue();
      aCompartment->lengthY = aSystem->getVariable("LENGTHY")->getValue();
      aCompartment->lengthZ = aSystem->getVariable("LENGTHZ")->getValue();
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
  //The center with reference to the immediate super system:
  aCompartment->centerPoint = aSuperCompartment->centerPoint;
  aCompartment->centerPoint.x += anOriX*aSuperCompartment->lengthX/2;
  aCompartment->centerPoint.y += anOriY*aSuperCompartment->lengthY/2;
  aCompartment->centerPoint.z += anOriZ*aSuperCompartment->lengthZ/2;
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


Point SpatiocyteStepper::coord2point(unsigned int aCoord)
{
  unsigned int aGlobalCol;
  unsigned int aGlobalLayer;
  unsigned int aGlobalRow;
  coord2global(aCoord, &aGlobalRow, &aGlobalLayer, &aGlobalCol);
  //the center point of a voxel 
  Point aPoint;
  //Specific for HCP lattice:{
  aPoint.y = (aGlobalCol%2)*theHCPk+theHCPl*aGlobalLayer;
  aPoint.z = aGlobalRow*2*theNormalizedVoxelRadius+
    ((aGlobalLayer+aGlobalCol)%2)*theNormalizedVoxelRadius;
  aPoint.x = aGlobalCol*theHCPh;
  //}
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
      if(getModel()->getRootSystem()->getVariable("SURFACEZ")->getValue() ==
         UNIPERIODIC)
        { 
          replaceUniVoxel(aSrcVoxel, aDestVoxel);
        }
      else if(getModel()->getRootSystem()->getVariable("SURFACEZ")->getValue()
              == PERIODIC)
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
      if(getModel()->getRootSystem()->getVariable("SURFACEY")->getValue() ==
         UNIPERIODIC)
        {
          replaceUniVoxel(aSrcVoxel, aDestVoxel);
        }
      else if(getModel()->getRootSystem()->getVariable("SURFACEY")->getValue()
              == PERIODIC)
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
      if(getModel()->getRootSystem()->getVariable("SURFACEX")->getValue() ==
         UNIPERIODIC)
        {
          replaceUniVoxel(aSrcVoxel, aDestVoxel);
        }
      else if(getModel()->getRootSystem()->getVariable("SURFACEX")->getValue()
              == PERIODIC)
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
  for(vector<Voxel>::iterator i(theLattice.begin()); i != theLattice.end(); ++i)
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
  for(vector<Compartment*>::iterator i(theCompartments.begin());
      i != theCompartments.end(); ++i)
    {
      if((*i)->isSurface)
        {
          removePeriodicEdgeVoxels(*i);
          for(vector<unsigned int>::iterator j((*i)->coords.begin());
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
      vector<unsigned int> coords;
      for(vector<unsigned int>::iterator j(aCompartment->coords.begin());
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

void SpatiocyteStepper::optimizeSurfaceVoxel(Voxel* aVoxel,
                                             Compartment* aCompartment)
{
  unsigned short surfaceID(aCompartment->vacantID);
  aVoxel->surfaceVoxels = new vector<vector<Voxel*> >;
  aVoxel->surfaceVoxels->resize(4);
  vector<Voxel*>& immediateSurface((*aVoxel->surfaceVoxels)[IMMEDIATE]);
  vector<Voxel*>& extendedSurface((*aVoxel->surfaceVoxels)[EXTENDED]);
  vector<Voxel*>& innerVolume((*aVoxel->surfaceVoxels)[INNER]);
  vector<Voxel*>& outerVolume((*aVoxel->surfaceVoxels)[OUTER]);
  vector<vector<Voxel*> > sharedVoxelsList;
  Voxel** forward(aVoxel->adjoiningVoxels);
  Voxel** reverse(forward+ADJOINING_VOXEL_SIZE);
  vector<Voxel*> adjoiningCopy;
  for(int k(0); k != ADJOINING_VOXEL_SIZE; ++k)
    {
      adjoiningCopy.push_back(forward[k]);
    }
  //Separate adjoining surface voxels and adjoining volume voxels.
  //Put the adjoining surface voxels at the beginning of the
  //adjoiningVoxels list while the volume voxels are put at the endl:
  for(vector<Voxel*>::iterator l(adjoiningCopy.begin());
      l != adjoiningCopy.end(); ++l)
    {
      if((*l)->id == surfaceID && (*l) != aVoxel)
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
                 find(adjoiningCopy.begin(), adjoiningCopy.end(),
                      extendedVoxel) == adjoiningCopy.end())
                {
                  vector<Voxel*>::iterator n(find(extendedSurface.begin(),
                        extendedSurface.end(), extendedVoxel));
                  if(n == extendedSurface.end())
                    {
                      extendedSurface.push_back(extendedVoxel);
                      //We require shared immediate voxel which
                      //connects the extended voxel with the source voxel 
                      //for polymerization. Create a new list of shared
                      //immediate voxel each time a new extended voxel is added:
                      vector<Voxel*> sharedVoxels;
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
  for(vector<vector<Voxel*> >::iterator i(sharedVoxelsList.begin());
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
  for(vector<Voxel>::iterator i(theLattice.begin()); i != theLattice.end(); ++i)
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
  if(aCompartment->isSurface)
    {
      if(isSurfaceVoxel(aVoxel, aCompartment))
        {
          aVoxel->id = aCompartment->vacantID;
          aCompartment->coords.push_back(aVoxel->coord-theStartCoord);
          return true;
        }
      return false;
    }
  else
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
          aVoxel->id = aCompartment->vacantID;
          aCompartment->coords.push_back(aVoxel->coord-theStartCoord);
          return true;
        }
      return false;
    }
}

bool SpatiocyteStepper::isSurfaceVoxel(Voxel* aVoxel, Compartment* aCompartment)
{
  if(isInsideCoord(aVoxel->coord, aCompartment, 2))
    {
      return false;
    }
  if(aCompartment->system->getSuperSystem()->isRootSystem())
    {
      for(unsigned int i(0); i != ADJOINING_VOXEL_SIZE; ++i)
        {
          if(aVoxel->adjoiningVoxels[i]->id == theNullID ||
             aVoxel->adjoiningVoxels[i] == aVoxel)
            {
              return true;
            }
        }
    }
  else
    {
      for(unsigned int i(0); i != ADJOINING_VOXEL_SIZE; ++i)
        {
          if(!isInsideCoord(aVoxel->adjoiningVoxels[i]->coord, aCompartment, 0))
            {
              return true;
            }
        }
    }
  return false;
}

bool SpatiocyteStepper::isInsideCoord(unsigned int aCoord,
                                      Compartment* aCompartment, double delta)
{
  Point aPoint(coord2point(aCoord));
  Point aCenterPoint(aCompartment->centerPoint);
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
  for(vector<Species*>::const_iterator i(aCompartment->species.begin());
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
          cout << "here" << endl;
          gsl_ran_choose(getRng(), populateVoxels, populationSize,
                     availableVoxels, availableVoxelSize, sizeof(unsigned int));
          //gsl_ran_choose arranges the position ascending, so we need
          //to shuffle the order of voxel positions:
          gsl_ran_shuffle(getRng(), populateVoxels, populationSize,
                          sizeof(unsigned int)); 
          for(vector<Species*>::const_iterator i(aCompartment->species.begin());
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
          for(vector<Species*>::const_iterator i(aCompartment->species.begin());
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
  for(vector<Species*>::const_iterator i(aCompartment->species.begin());
      i != aCompartment->species.end(); ++i)
    {
      (*i)->removeMolecules();
    }
}




