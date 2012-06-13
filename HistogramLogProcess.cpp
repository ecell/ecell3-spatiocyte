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

#include "HistogramLogProcess.hpp"

void HistogramLogProcess::initializeFifth()
{
  for(std::vector<Species*>::const_iterator i(theProcessSpecies.begin());
      i != theProcessSpecies.end(); ++i)
    {
      if((*i)->getDiffusionInterval() < theStepInterval)
        {
          theStepInterval = (*i)->getDiffusionInterval();
        }
    }
  if(LogInterval > 0)
    {
      theStepInterval = LogInterval;
    }
  else
    {
      LogInterval = theStepInterval;
    }
  theTime = LogInterval;
  thePriorityQueue->move(theQueueID);
}

void HistogramLogProcess::initializeLastOnce()
{
  theLogFile.open(FileName.c_str(), std::ios::trunc);
  theTotalIterations = Iterations;
  timePoints = (unsigned int)ceil(LogDuration/theStepInterval);
  theLogValues.resize(timePoints);
  for(unsigned int i(0); i != timePoints; ++i)
    {
      theLogValues[i].resize(Bins);
      for(unsigned int j(0); j != Bins; ++j)
        {
          theLogValues[i][j].resize(theProcessSpecies.size());
          for(unsigned int k(0); k != theProcessSpecies.size(); ++k)
            {
              theLogValues[i][j][k] = 0;
            }
        }
    }
}

void HistogramLogProcess::fire()
{
  if(theTime <= LogDuration)
    {
      logValues();
      ++timePointCnt;
    }
  if(theTime >= LogDuration && Iterations > 0)
    {
      theStepInterval = LogInterval;
      --Iterations;
      std::cout << "Iterations left:" << Iterations << " of " <<
        theTotalIterations << std::endl;
      if(Iterations > 0)
        {
          theSpatiocyteStepper->reset(Iterations);
          saveBackup();
          return;
        }
    }
  if(Iterations == 0)
    {
      saveFile();
      std::cout << "Done saving." << std::endl;
    }
  theTime += theStepInterval;
  thePriorityQueue->moveTop();
}


void HistogramLogProcess::saveFile()
{
  std::cout << "Saving data in: " << FileName.c_str() << std::endl;
  double aTime(LogInterval);
  for(unsigned int i(0); i != timePoints; ++i)
    {
      for(unsigned int j(0); j != Bins; ++j)
        {
          theLogFile << std::setprecision(15) << aTime << "," << j;
          for(unsigned int k(0); k != theProcessSpecies.size(); ++k)
            {
              theLogFile << "," << std::setprecision(15) <<
                theLogValues[i][j][k]/theTotalIterations;
            }
          theLogFile << std::endl;
        }
      aTime += LogInterval;
    }
  theLogFile.close();
  theStepInterval = libecs::INF;
}

void HistogramLogProcess::saveBackup()
{
  if(SaveCounts > 0 && 
     Iterations%(int)rint(theTotalIterations/SaveCounts) == 0)
    {
      std::string aFileName(FileName.c_str());
      aFileName = aFileName + ".back";
      std::cout << "Saving backup data in: " << aFileName 
        << std::endl;
      std::ofstream aFile;
      aFile.open(aFileName.c_str(), std::ios::trunc);
      double aTime(LogInterval);
      int completedIterations(theTotalIterations-Iterations);
      for(unsigned int i(0); i != timePoints; ++i)
        {
          for(unsigned int j(0); j != Bins; ++j)
            {
             aFile << std::setprecision(15) << aTime << "," << j;
              for(unsigned int k(0); k != theProcessSpecies.size(); ++k)
                {
                  aFile << "," << std::setprecision(15) <<
                    theLogValues[i][j][k]/completedIterations;
                }
              aFile << std::endl;
            }
          aTime += LogInterval;
        }
      aFile.close();
    }
}

void HistogramLogProcess::logValues()
{
 //std::cout << "timePoint:" << timePointCnt <<  " curr:" << theSpatiocyteStepper->getCurrentTime() << std::endl;
  for(unsigned int i(0); i != theProcessSpecies.size(); ++i)
    {
      Species* aSpecies(theProcessSpecies[i]);
      std::cout << getIDString(aSpecies) << " " << aSpecies->size() << std::endl;
      for(unsigned int j(0); j != aSpecies->size(); ++j)
        {
          unsigned int bin;
          if(isInside(bin, aSpecies->getPoint(j)))
            {
              ++theLogValues[timePointCnt][bin][i];
            }
        }
    }
}

void HistogramLogProcess::initializeFourth()
{
  theComp = theSpatiocyteStepper->system2Comp(getSuperSystem());
  C = theComp->centerPoint;
  C.x += OriginX*theComp->lengthX/2;
  C.y += OriginY*theComp->lengthY/2;
  C.z += OriginZ*theComp->lengthZ/2;
  VoxelDiameter = theSpatiocyteStepper->getVoxelRadius()*2;
  Length /= VoxelDiameter;
  Radius /= VoxelDiameter;
  binInterval = Length/Bins;
  initializeVectors();
}

void HistogramLogProcess::initializeVectors()
{
  E.x = -Length/2;
  E.y = 0;
  E.z = 0;
  //Rotated Minus end
  theSpatiocyteStepper->rotateX(RotateX, &E, -1);
  theSpatiocyteStepper->rotateY(RotateY, &E, -1);
  theSpatiocyteStepper->rotateZ(RotateZ, &E, -1);
  E.x += C.x;
  E.y += C.y;
  E.z += C.z;
  //Direction vector from the East to center
  D.x = C.x-E.x;
  D.y = C.y-E.y;
  D.z = C.z-E.z;
  //Make D a unit vector
  double NormT(sqrt(D.x*D.x+D.y*D.y+D.z*D.z));
  D.x /= NormT;
  D.y /= NormT;
  D.z /= NormT;
}


bool HistogramLogProcess::isInside(unsigned int& bin, Point N)
{
  double t((D.x*N.x+D.y*N.y+D.z*N.z-D.x*E.x-D.y*E.y-D.z*E.z)/
           (D.x*D.x +D.y*D.y+D.z*D.z));
  if(t < 0)
    {
      return false;
    }
  bin = (unsigned int)floor(t/binInterval);
  if(bin >= Bins)
    {
      return false;
    }
  double dist(sqrt(pow(-N.x+E.x+D.x*t, 2)+pow(-N.y+E.y+D.y*t, 2)+
                   pow(-N.z+E.z+D.z*t, 2)));
  if(dist > Radius)
    {
      return false;
    }
  return true;
}

LIBECS_DM_INIT(HistogramLogProcess, Process); 
