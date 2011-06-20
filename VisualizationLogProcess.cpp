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

#include "VisualizationLogProcess.hpp"

LIBECS_DM_INIT(VisualizationLogProcess, Process); 

void VisualizationLogProcess::initializeLog()
{
  unsigned int aThreadSize(1);
  unsigned int aLatticeType(theSpatiocyteStepper->getLatticeType());
  theLogFile.write((char*)(&aThreadSize), sizeof(aThreadSize));
  theLogFile.write((char*)(&aLatticeType), sizeof(aLatticeType));
  theLogFile.write((char*)(&theMeanCount), sizeof(theMeanCount));
  unsigned int aStartCoord(theSpatiocyteStepper->getStartCoord());
  theLogFile.write((char*)(&aStartCoord), sizeof(aStartCoord));
  unsigned int aRowSize(theSpatiocyteStepper->getRowSize());
  theLogFile.write((char*)(&aRowSize), sizeof(aRowSize));
  unsigned int aLayerSize(theSpatiocyteStepper->getLayerSize());
  theLogFile.write((char*)(&aLayerSize), sizeof(aLayerSize));
  unsigned int aColSize(theSpatiocyteStepper->getColSize());
  theLogFile.write((char*)(&aColSize), sizeof(aColSize));
  Point aCenterPoint(theSpatiocyteStepper->getCenterPoint());
  double aRealRowSize(aCenterPoint.z*2);
  theLogFile.write((char*)(&aRealRowSize), sizeof(aRealRowSize));
  double aRealLayerSize(aCenterPoint.y*2);
  theLogFile.write((char*)(&aRealLayerSize), sizeof(aRealLayerSize));
  double aRealColSize(aCenterPoint.x*2);
  theLogFile.write((char*)(&aRealColSize), sizeof(aRealColSize));
  unsigned int aSpeciesSize(theProcessSpecies.size());
  theLogFile.write((char*)(&aSpeciesSize), sizeof(aSpeciesSize));
  unsigned int aPolymerSize(thePolymerSpecies.size());
  theLogFile.write((char*)(&aPolymerSize), sizeof(aPolymerSize));
  unsigned int aReservedSize(0);
  theLogFile.write((char*)(&aReservedSize), sizeof(aReservedSize));
  //theLogMarker is a constant throughout the simulation:
  theLogFile.write((char*)(&theLogMarker), sizeof(theLogMarker));
  double aVoxelRadius(theSpatiocyteStepper->getNormalizedVoxelRadius());
  theLogFile.write((char*)(&aVoxelRadius), sizeof(aVoxelRadius));
  for(unsigned int i(0); i != theProcessSpecies.size(); ++i)
    {
      unsigned int aStringSize(
       theProcessSpecies[i]->getVariable()->getFullID().asString().size());
      theLogFile.write((char*)(&aStringSize), sizeof(aStringSize));
      theLogFile.write(
       theProcessSpecies[i]->getVariable()->getFullID().asString().c_str(),
       aStringSize);
    }
  for(unsigned int i(0); i!=thePolymerSpecies.size(); ++i)
    {
      unsigned int aPolymerIndex(thePolymerIndex[i]);
      theLogFile.write((char*) (&aPolymerIndex), sizeof(aPolymerIndex));
    }
  //a, b, c are used for multithreaded simulation which is
  //not implemented yet.
  //The visualizer assumes multithreaded data, so we need to initialize the 
  //logger with the expected values for a single threaded simulation.
  unsigned int a(aStartCoord-1);
  unsigned int c(a+aRowSize*aLayerSize*aColSize);
  unsigned int b(c-aRowSize*aLayerSize*0);
  theLogFile.write((char*)(&a), sizeof(a));
  theLogFile.write((char*)(&b), sizeof(b));
  theLogFile.write((char*)(&c), sizeof(c));
  theLogFile.flush();
}

void VisualizationLogProcess::logMolecules(int anIndex)
{
  Species* aSpecies(theProcessSpecies[anIndex]);
  //No need to log lipid or vacant molecules since the size is 0:
  if(aSpecies->getIsLipid() || aSpecies->getIsVacant())
    {
      return;
    }
  theLogFile.write((char*)(&anIndex), sizeof(anIndex));
  //The species molecule size:
  int aSize(aSpecies->size());
  theLogFile.write((char*)(&aSize), sizeof(aSize)); 
  for(int i(0); i != aSize; ++i)
    {
      unsigned int aCoord(aSpecies->getCoord(i));
      theLogFile.write((char*)(&aCoord), sizeof(aCoord));
    }
}  

void VisualizationLogProcess::logPolymers(int anIndex)
{
  Species* aSpecies(thePolymerSpecies[anIndex]);
  theLogFile.write((char*)(&anIndex), sizeof(anIndex));
  //The species molecule size:
  int aSize(aSpecies->size());
  theLogFile.write((char*)(&aSize), sizeof(aSize)); 
  for(int i(0); i != aSize; ++i)
    {
      Point aPoint(aSpecies->getPoint(i));
      theLogFile.write((char*)(&aPoint), sizeof(aPoint));
    }
}  

void VisualizationLogProcess::logSourceMolecules(int anIndex)
{
  Species* aSpecies(thePolymerSpecies[anIndex]);
  int aSourceIndex(theProcessSpecies.size()+anIndex);
  theLogFile.write((char*)(&aSourceIndex), sizeof(aSourceIndex));
  const std::vector<unsigned int> aCoords(aSpecies->getSourceCoords());
  int aSize(aCoords.size());
  theLogFile.write((char*)(&aSize), sizeof(aSize)); 
  for(unsigned int i(0); i != aCoords.size(); ++i)
    {
      unsigned int aCoord(aCoords[i]);
      theLogFile.write((char*)(&aCoord), sizeof(aCoord));
    }
}  

void VisualizationLogProcess::logTargetMolecules(int anIndex)
{
  Species* aSpecies(thePolymerSpecies[anIndex]);
  int aTargetIndex(theProcessSpecies.size()+thePolymerSpecies.size()+anIndex);
  theLogFile.write((char*)(&aTargetIndex), sizeof(aTargetIndex));
  const std::vector<unsigned int> aCoords(aSpecies->getTargetCoords());
  int aSize(aCoords.size());
  theLogFile.write((char*)(&aSize), sizeof(aSize)); 
  for(unsigned int i(0); i != aCoords.size(); ++i)
    {
      unsigned int aCoord(aCoords[i]);
      theLogFile.write((char*)(&aCoord), sizeof(aCoord));
    }
}  

void VisualizationLogProcess::logSharedMolecules(int anIndex)
{
  Species* aSpecies(thePolymerSpecies[anIndex]);
  int aSharedIndex(theProcessSpecies.size()+thePolymerSpecies.size()*2+anIndex);
  theLogFile.write((char*)(&aSharedIndex), sizeof(aSharedIndex));
  const std::vector<unsigned int> aCoords(aSpecies->getSharedCoords());
  int aSize(aCoords.size());
  theLogFile.write((char*)(&aSize), sizeof(aSize)); 
  for(unsigned int i(0); i != aCoords.size(); ++i)
    {
      unsigned int aCoord(aCoords[i]);
      theLogFile.write((char*)(&aCoord), sizeof(aCoord));
    }
}  

void VisualizationLogProcess::logSpecies()
{
  int aDataSize(0);
  std::streampos aStartPos(theLogFile.tellp());
  // write the next size (create a temporary space for it) 
  theLogFile.write((char*)(&aDataSize), sizeof(aDataSize));
  double aCurrentTime(theSpatiocyteStepper->getCurrentTime());
  theLogFile.write((char*)(&aCurrentTime), sizeof(aCurrentTime));
  for(unsigned int i(0); i != theProcessSpecies.size(); ++i)
    {
      logMolecules(i);
    }
  for(unsigned int i(0); i != thePolymerSpecies.size(); ++i)
    {
      logSourceMolecules(i);
    }
  for(unsigned int i(0); i != thePolymerSpecies.size(); ++i)
    {
      logTargetMolecules(i);
    }
  for(unsigned int i(0); i != thePolymerSpecies.size(); ++i)
    {
      logSharedMolecules(i);
    }
  /*
  for(unsigned int i(0); i!=theReservedSpecies.size(); ++i)
    {
      logReservedMolecules(i);
    }
    */
  //theLogMarker is a constant throughout the simulation:
  theLogFile.write((char*)(&theLogMarker), sizeof(theLogMarker));
  for(unsigned int i(0); i != thePolymerSpecies.size(); ++i)
    {
      logPolymers(i);
    }
  //theLogMarker is a constant throughout the simulation:
  theLogFile.write((char*)(&theLogMarker), sizeof(theLogMarker));
  aDataSize = (theLogFile.tellp()-aStartPos)-static_cast<std::streampos>(sizeof(aDataSize));
  int aPrevDataSize(theLogFile.tellp()-theStepStartPos+static_cast<std::streampos>(sizeof(int))*2);
  theStepStartPos = theLogFile.tellp();
  // write the prev size at the end of this step
  theLogFile.write((char*)(&aPrevDataSize), sizeof(aPrevDataSize));
  std::streampos aCurrentPos(theLogFile.tellp());
  theLogFile.seekp(aStartPos);
  // write the next size at the beginning of this step
  theLogFile.write((char*) (&aDataSize), sizeof(aDataSize));
  theLogFile.seekp(aCurrentPos);
}

void VisualizationLogProcess::logSurfaceVoxels()
{
  double aCurrentTime(theSpatiocyteStepper->getCurrentTime());
  theStepStartPos = theLogFile.tellp();
  // write prev pos = 0;
  int aDataSize(sizeof(int)*2);
  theLogFile.write((char*)(&aDataSize), sizeof(aDataSize));
  std::streampos aStartPos(theLogFile.tellp());
  // write the next size (create a temporary space for it) 
  theLogFile.write((char*)(&aDataSize), sizeof(aDataSize));
  theLogFile.write((char*)(&aCurrentTime), sizeof(aCurrentTime));
  for(unsigned int i(0); i != theProcessSpecies.size(); ++i)
    {
      if(theProcessSpecies[i]->getIsLipid())
        {
          Species* aLipid(theProcessSpecies[i]);
          //The species index in the process:
          theLogFile.write((char*)(&i), sizeof(i));
          const Comp* aSurface(aLipid->getComp());
          //The species molecule size:
          unsigned int aSize(aSurface->coords.size());
          theLogFile.write((char*)(&aSize), sizeof(aSize)); 
          unsigned int aStartCoord(theSpatiocyteStepper->getStartCoord());
          for(std::vector<unsigned int>::const_iterator j(
               aSurface->coords.begin()); j != aSurface->coords.end(); ++j)
            {
              unsigned int aCoord((*j)+aStartCoord);
              theLogFile.write((char*)(&aCoord), sizeof(aCoord));
            }  
        }
    }
  //theLogMarker is a constant throughout the simulation:
  theLogFile.write((char*)(&theLogMarker), sizeof(theLogMarker));
  theLogFile.write((char*)(&theLogMarker), sizeof(theLogMarker));
  aDataSize = (theLogFile.tellp()-aStartPos)-static_cast<std::streampos>(sizeof(aDataSize)); 
  int aPrevDataSize(theLogFile.tellp()-theStepStartPos+static_cast<std::streampos>(sizeof(int))*2);
  theStepStartPos = theLogFile.tellp();
  // write the prev size at the end of this step
  theLogFile.write((char*)(&aPrevDataSize), sizeof(aPrevDataSize));
  std::streampos aCurrentPos(theLogFile.tellp());
  theLogFile.seekp(aStartPos);
  // write the next size at the beginning of this step
  theLogFile.write((char*)(&aDataSize), sizeof(aDataSize));
  theLogFile.seekp(aCurrentPos);
}
