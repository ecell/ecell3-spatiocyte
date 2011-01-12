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

#include "FluorescentImagingProcess.hpp"

LIBECS_DM_INIT(FluorescentImagingProcess, Process); 

void FluorescentImagingProcess::incSpeciesLatticeCount()
{
  for(unsigned int i(0); i != thePositiveSpecies.size(); ++i)
    {
      Species* aSpecies(thePositiveSpecies[i]);
      unsigned int aMoleculeSize(aSpecies->size());
      for(unsigned int j(0); j != aMoleculeSize; ++j)
        { 
          Voxel* aMolecule(aSpecies->getMolecule(j));
          ++theLattice[theProcessSpeciesIndices[i][0]][
            aMolecule->coord-theStartCoord];
          for(unsigned int k(1); k != theProcessSpeciesIndices[i].size(); ++k)
            {
              Voxel* anAdjoin(aMolecule->adjoiningVoxels[k-1]);
              ++theLattice[theProcessSpeciesIndices[i][k]][
                anAdjoin->coord-theStartCoord];
            }
        }
    }
}

void FluorescentImagingProcess::logFluorescentSpecies()
{
  std::vector<int> coordList;
  for(unsigned int i(0); i != theLatticeSize; ++i)
    { 
      for(unsigned int j(0); j != theProcessSpecies.size(); ++j)
        {
          if(theLattice[j][i])
            {
              coordList.push_back(i);
              break;
            }
        }
    }
  int aDataSize(0);
  std::streampos aStartPos(theLogFile.tellp());
  // write the next size (create a temporary space for it) 
  theLogFile.write((char*)(&aDataSize), sizeof(aDataSize));
  double aCurrentTime(theSpatiocyteStepper->getCurrentTime());
  theLogFile.write((char*)(&aCurrentTime), sizeof(aCurrentTime));
  unsigned int coordSize(coordList.size());
  theLogFile.write((char*)(&coordSize), sizeof(coordSize));
  for(unsigned int i(0); i != coordSize; ++i)
    {
      unsigned int aCoord(coordList[i]+theStartCoord);
      theLogFile.write((char*)(&aCoord), sizeof(aCoord));
    }
  for(unsigned int i(0); i != theProcessSpecies.size(); ++i)
    {
      for(unsigned int j(0); j != coordSize; ++j)
        {
          unsigned int frequency(theLattice[i][coordList[j]]);
          theLogFile.write((char*)(&frequency), sizeof(frequency));
        }
    }
  aDataSize = (theLogFile.tellp()-aStartPos)-sizeof(aDataSize); 
  int aPrevDataSize(theLogFile.tellp()-theStepStartPos+sizeof(int)*2);
  theStepStartPos = theLogFile.tellp();
  // write the prev size at the end of this step
  theLogFile.write((char*)(&aPrevDataSize), sizeof(aPrevDataSize));
  std::streampos aCurrentPos(theLogFile.tellp());
  theLogFile.seekp(aStartPos);
  // write the next size at the beginning of this step
  theLogFile.write((char*) (&aDataSize), sizeof(aDataSize));
  theLogFile.seekp(aCurrentPos);
}

