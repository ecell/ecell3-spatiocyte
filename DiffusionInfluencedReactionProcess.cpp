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

#include "DiffusionInfluencedReactionProcess.hpp"
#include "SpatiocyteSpecies.hpp"

LIBECS_DM_INIT(DiffusionInfluencedReactionProcess, Process);


void DiffusionInfluencedReactionProcess::initializeSecond()
{
  ReactionProcess::initializeSecond();
  A->setDiffusionInfluencedReactantPair(B); 
  B->setDiffusionInfluencedReactantPair(A); 
  r_v = theSpatiocyteStepper->getVoxelRadius();
  D_A = A->getDiffusionCoefficient();
  D_B = B->getDiffusionCoefficient();
  calculateReactionProbability();
  if(A->getIsDiffusing())
    {
      A->setDiffusionInfluencedReaction(this, B->getID(), p); 
    }
  if(B->getIsDiffusing())
    {
      B->setDiffusionInfluencedReaction(this, A->getID(), p); 
    }
}

void DiffusionInfluencedReactionProcess::initializeThird()
{
}

//Do the reaction A + B -> C + D. So that A <- C and B <- D.
//We need to consider that the source molecule can be either A or B.
//If A and C belong to the same Comp, A <- C.
//Otherwise, find a vacant adjoining voxel of A, X which is the same Comp
//as C and X <- C.
//Similarly, if B and D belong to the same Comp, B <- D.
//Otherwise, find a vacant adjoining voxel of C, Y which is the same Comp
//as D and Y <- D.
bool DiffusionInfluencedReactionProcess::react(Voxel* moleculeB, Voxel** target)
{
  Voxel* moleculeA(*target);
  //moleculeA is the source molecule. It will be soft-removed (id kept intact)
  //by the calling Species if this method returns true.
  //moleculeB is the target molecule,  it will also be soft-removed by the
  //calling Species.
  //First let us make sure moleculeA and moleculeB belong to the
  //correct species.
  if(moleculeA->id != A->getID())
    {
      Voxel* tempA(moleculeA);
      moleculeA = moleculeB;
      moleculeB = tempA;
    }
  //nonHD_A + nonHD_B -> nonHD_C + HD_D:
  //nonHD_A + nonHD_B -> HD_C + nonHD_D:
  if((variableC && D) || (C && variableD))
    {
      Variable* HD_p(variableC);
      Species* nonHD_p(D);
      if(variableD)
        {
          HD_p = variableD;
          nonHD_p = C;
        }
      Voxel* moleculeP;
      if(A->getVacantID() == nonHD_p->getVacantID() ||
         A->getID() == nonHD_p->getVacantID())
        {
          moleculeP = moleculeA;
          //Hard remove the B molecule, since nonHD_p is in a different Comp:
          B->getVacantSpecies()->addMolecule(moleculeB);
        }
      else if(B->getVacantID() == nonHD_p->getVacantID() ||
              B->getID() == nonHD_p->getVacantID())
        {
          moleculeP = moleculeB;
          //Hard remove the A molecule, since nonHD_p is in a different Comp:
          A->getVacantSpecies()->addMolecule(moleculeA);
        }
      else
        { 
          moleculeP = nonHD_p->getRandomAdjoiningVoxel(moleculeA);
          //Only proceed if we can find an adjoining vacant voxel
          //of A which can be occupied by C:
          if(moleculeP == NULL)
            {
              moleculeP = nonHD_p->getRandomAdjoiningVoxel(moleculeB);
              if(moleculeP == NULL)
                {
                  return false;
                }
            }
          //Hard remove the A molecule, since nonHD_p is in a different Comp:
          A->getVacantSpecies()->addMolecule(moleculeA);
          //Hard remove the B molecule, since nonHD_p is in a different Comp:
          B->getVacantSpecies()->addMolecule(moleculeB);
        }
      HD_p->addValue(1);
      nonHD_p->addMolecule(moleculeP);
      return true;
    }
  //nonHD_A + nonHD_B -> HD_C:
  else if(variableC && !D && !variableD)
    {

      //Hard remove the A molecule, since nonHD_p is in a different Comp:
      A->getVacantSpecies()->addMolecule(moleculeA);
      //Hard remove the B molecule, since nonHD_p is in a different Comp:
      B->getVacantSpecies()->addMolecule(moleculeB);
      variableC->addValue(1);
      return true;
    }

  //If the product C is not in the same Comp as A,
  //we need to find a vacant adjoining voxel of A that belongs
  //to the Comp of C:
  Voxel* moleculeC;
  if(A->getVacantID() == C->getVacantID() || A->getID() == C->getVacantID())
    {
      moleculeC = moleculeA;
    }
  else
    {
      moleculeC = C->getRandomAdjoiningVoxel(moleculeA);
      //Only proceed if we can find an adjoining vacant voxel
      //of A which can be occupied by C:
      if(moleculeC == NULL)
        {
          return false;
        }
    }
  //If it has two products:
  if(D != NULL)
    { 
      //If the product D is not in the same Comp as B,
      //we need to find a vacant adjoining voxel of C that belongs
      //to the Comp of D:
      Voxel* moleculeD;
      if(B->getVacantID() == D->getVacantID() || B->getID() == D->getVacantID())
        {
          moleculeD = moleculeB;
        }
      else
        {
          moleculeD = D->getRandomAdjoiningVoxel(moleculeC, moleculeC);
          //Only proceed if we can find an adjoining vacant voxel
          //of A which can be occupied by C:
          if(moleculeD == NULL)
            {
              return false;
            }
          B->getVacantSpecies()->addMolecule(moleculeB);
        }
      D->addMolecule(moleculeD);
    }
  else
    {
      //Hard remove the B molecule since this is a single product reaction:
      B->getVacantSpecies()->addMolecule(moleculeB);
    }
  //Hard remove the A molecule, in case C is in a different Comp:
  A->getVacantSpecies()->addMolecule(moleculeA);
  C->addMolecule(moleculeC);
}


  

void DiffusionInfluencedReactionProcess::finalizeReaction()
{
  //The number of molecules may have changed for both reactant and product
  //species. We need to update SpatiocyteNextReactionProcesses which are
  //dependent on these species:
  for(std::vector<ReactionProcess*>::const_iterator 
      i(theInterruptingProcesses.begin());
      i!=theInterruptingProcesses.end(); ++i)
    {
      (*i)->substrateValueChanged(theSpatiocyteStepper->getCurrentTime());
    }
}

void DiffusionInfluencedReactionProcess::calculateReactionProbability()
{
  //Refer to the paper for the description of the variables used in this
  //method.
  if(A->getIsVolume() && B->getIsVolume())
    {
      if(A != B)
        {
          if(p == -1)
            {
              p = k/(6*sqrt(2)*(D_A+D_B)*r_v);
            }
          else
            {
              k = p*(6*sqrt(2)*(D_A+D_B)*r_v);
            }
        }
      else
        {
          if(p == -1)
            {
              p = k/(6*sqrt(2)*D_A*r_v);
            }
          else
            {
              k = p*(6*sqrt(2)*D_A*r_v);
            }
        }
    }
  else if(!A->getIsVolume() && !B->getIsVolume())
    {
      //Inter-surface Comp reaction.
      //For surface edge absorbing reactions:
      if(A->getComp() != B->getComp())
        {
          k = p;
        }
      else if(A != B)
        {
          if(p == -1)
            {
              p = pow(2*sqrt(2)+4*sqrt(3)+3*sqrt(6)+sqrt(22), 2)*k/
                (72*(6*sqrt(2)+4*sqrt(3)+3*sqrt(6))*(D_A+D_B));
            }
          else
            {
              k = p*(72*(6*sqrt(2)+4*sqrt(3)+3*sqrt(6))*(D_A+D_B))/
                pow(2*sqrt(2)+4*sqrt(3)+3*sqrt(6)+sqrt(22), 2);
            }
        }
      else
        {
          if(p == -1)
            {
              p = pow(2*sqrt(2)+4*sqrt(3)+3*sqrt(6)+sqrt(22), 2)*k/
                (72*(6*sqrt(2)+4*sqrt(3)+3*sqrt(6))*(D_A));
            }
          else
            {
              k = p*(72*(6*sqrt(2)+4*sqrt(3)+3*sqrt(6))*(D_A))/
                pow(2*sqrt(2)+4*sqrt(3)+3*sqrt(6)+sqrt(22), 2);
            }
        }
    }
  else if(A->getIsVolume() && B->getIsLipid())
    {
      if(p == -1)
        {
          p = 24*k*r_v/((6+3*sqrt(3)+2*sqrt(6))*D_A);
        }
      else
        {
          k = p*((6+3*sqrt(3)+2*sqrt(6))*D_A)/(24*r_v);
        }
    }
  else if(A->getIsLipid() && B->getIsVolume())
    {
      if(p == -1)
        {
          p = 24*k*r_v/((6+3*sqrt(3)+2*sqrt(6))*D_B);
        }
      else
        {
          k = p*((6+3*sqrt(3)+2*sqrt(6))*D_B)/(24*r_v);
        }
    }
  else if(A->getIsVolume() && !B->getIsVolume())
    {
      if(p == -1)
        {
          p = sqrt(2)*k/(3*D_A*r_v);
        }
      else
        {
          k = p*(3*D_A*r_v)/sqrt(2);
        }
    }
  else if(!A->getIsVolume() && B->getIsVolume())
    {
      if(p == -1)
        {
          p = sqrt(2)*k/(3*D_B*r_v);
        }
      else
        {
          k = p*(3*D_B*r_v)/sqrt(2);
        }
    }
  else
    {
      THROW_EXCEPTION(ValueError, 
                      String(getPropertyInterface().getClassName()) + 
                      "[" + getFullID().asString() + 
                      "]: Error in type of second order reaction.");
    }
}

void DiffusionInfluencedReactionProcess::printParameters()
{
  String aProcess(String(getPropertyInterface().getClassName()) + 
                                      "[" + getFullID().asString() + "]");
  std::cout << aProcess << std::endl;
  std::cout << "  " << getIDString(A) << " + " <<  getIDString(B) << " -> ";
  if(C)
    {
      std::cout << getIDString(C);
    }
  else
    {
      std::cout << getIDString(variableC);
    }
  if(D)
    {
      std::cout << " + " << getIDString(D);
    }
  else if(variableD)
    {
      std::cout << " + " << getIDString(variableD);
    }
  std::cout << ": k=" << k << ", p=" << p << 
    ", p_A=" << A->getReactionProbability(B->getID()) <<
    ", p_B=" << B->getReactionProbability(A->getID()) << std::endl; 
}


    
      


