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

#include "MicroscopyTrackingProcess.hpp"

LIBECS_DM_INIT(MicroscopyTrackingProcess, Process); 

void MicroscopyTrackingProcess::incSpeciesLatticeCount()
{
  for(unsigned int i(0); i != thePositiveSpecies.size(); ++i)
    {
      Species* aSpecies(thePositiveSpecies[i]);
      unsigned int aMoleculeSize(aSpecies->size());
      for(unsigned int j(0); j != aMoleculeSize; ++j)
        { 
          Voxel* aMolecule(aSpecies->getMolecule(j));
          ++theLattice[theLatticeSpeciesIndices[i][0]][
            aMolecule->coord-theStartCoord];
          for(unsigned int k(1); k != theLatticeSpeciesIndices[i].size(); ++k)
            {
              Voxel* anAdjoin(aMolecule->adjoiningVoxels[k-1]);
              ++theLattice[theLatticeSpeciesIndices[i][k]][
                anAdjoin->coord-theStartCoord];
            }
        }
    }
}

void MicroscopyTrackingProcess::logFluorescentSpecies()
{
  std::vector<int> coordList;
  for(unsigned int i(0); i != theLatticeSize; ++i)
    { 
      for(unsigned int j(0); j != theLatticeSpecies.size(); ++j)
        {
          if(theLattice[j][i])
            {
              coordList.push_back(i);
              break;
            }
        }
    }
  double aCurrentTime(theSpatiocyteStepper->getCurrentTime());
  theLogFile.write((char*)(&aCurrentTime), sizeof(aCurrentTime));
  unsigned int coordSize(coordList.size());
  theLogFile.write((char*)(&coordSize), sizeof(coordSize));
  for(unsigned int i(0); i != coordSize; ++i)
    {
      unsigned int aCoord(coordList[i]+theStartCoord);
      theLogFile.write((char*)(&aCoord), sizeof(aCoord));
    }
  for(unsigned int i(0); i != theLatticeSpecies.size(); ++i)
    {
      for(unsigned int j(0); j != coordSize; ++j)
        {
          unsigned int frequency(theLattice[i][coordList[j]]);
          theLogFile.write((char*)(&frequency), sizeof(frequency));
        }
    }
}

