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

void DiffusionInfluencedReactionProcess::checkSubstrates()
{
  //HD_A or HD_B:
	if(variableA)
    {
      THROW_EXCEPTION(ValueError, String(
        getPropertyInterface().getClassName()) + " [" + getFullID().asString() +
		    "]: A DiffusionInfluencedReactionProcess cannot have a HD " +
        "substrate species: " + getIDString(variableA));
    }
  if(variableB)
    {
      THROW_EXCEPTION(ValueError, String(
        getPropertyInterface().getClassName()) + " [" + getFullID().asString() +
		    "]: A DiffusionInfluencedReactionProcess cannot have a HD " +
        "substrate species: " + getIDString(variableB));
    }
}

void DiffusionInfluencedReactionProcess::initializeSecond()
{
  ReactionProcess::initializeSecond(); 
  A->setCollision(Collision);
  B->setCollision(Collision);
}

void DiffusionInfluencedReactionProcess::initializeThird()
{
  ReactionProcess::initializeThird();
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

//Do the reaction A + B -> C + D. So that A <- C and B <- D.
//We need to consider that the source molecule can be either A or B.
//If A and C belong to the same Comp, A <- C.
//Otherwise, find a vacant adjoining voxel of A, X which is the same Comp
//as C and X <- C.
//Similarly, if B and D belong to the same Comp, B <- D.
//Otherwise, find a vacant adjoining voxel of C, Y which is the same Comp
//as D and Y <- D.
bool DiffusionInfluencedReactionProcess::react(unsigned int coordA,
                                               unsigned int coordB)
{
  //First let us make sure moleculeA and moleculeB belong to the
  //correct species.
  if((*theLattice)[coordA].id != A->getID())
    {
      unsigned int tempA(coordA);
      coordA = coordB;
      coordB = tempA;
    }
  Voxel& moleculeA((*theLattice)[coordA]);
  Voxel& moleculeB((*theLattice)[coordB]);

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
      unsigned int coordP;
      if(A->getVacantID() == nonHD_p->getVacantID() ||
         A->getID() == nonHD_p->getVacantID())
        {
          coordP = coordA;
          //Hard remove the B molecule, since nonHD_p is in a different Comp:
          moleculeB.id = B->getVacantID();
        }
      else if(B->getVacantID() == nonHD_p->getVacantID() ||
              B->getID() == nonHD_p->getVacantID())
        {
          coordP = coordB;
          //Hard remove the A molecule, since nonHD_p is in a different Comp:
          moleculeA.id = A->getVacantID();
        }
      else
        { 
          coordP = nonHD_p->getRandomAdjoiningCoord(coordA, SearchVacant);
          //Only proceed if we can find an adjoining vacant voxel
          //of A which can be occupied by C:
          if(coordP == theNullCoord)
            {
              coordP = nonHD_p->getRandomAdjoiningCoord(coordB, SearchVacant);
              if(coordP == theNullCoord)
                {
                  return false;
                }
            }
          //Hard remove the A molecule, since nonHD_p is in a different Comp:
          moleculeA.id = A->getVacantID();
          //Hard remove the B molecule, since nonHD_p is in a different Comp:
          moleculeB.id = B->getVacantID();
        }
      HD_p->addValue(1);
      nonHD_p->addCoord(coordP);
      return true;
    }
  //nonHD_A + nonHD_B -> HD_C:
  else if(variableC && !D && !variableD)
    {

      //Hard remove the A molecule, since nonHD_p is in a different Comp:
      moleculeA.id = A->getVacantID();
      //Hard remove the B molecule, since nonHD_p is in a different Comp:
      moleculeB.id = B->getVacantID();
      variableC->addValue(1);
      return true;
    }
  unsigned int coordC;
  unsigned int coordD;
  if(A->getVacantID() == C->getVacantID() || A->getID() == C->getVacantID())
    {
      coordC = coordA;
      if(D)
        {
          if(B->getVacantID() == D->getVacantID() ||
             B->getID() == D->getVacantID())
            {
              coordD = coordB;
            }
          else
            {
              coordD = D->getRandomAdjoiningCoord(coordC, coordC, SearchVacant);
              if(coordD == theNullCoord)
                {
                  return false;
                }
              moleculeB.id = B->getVacantID();
            }
          D->addCoord(coordD);
        }
      else
        {
          //Hard remove the B molecule since it is not used:
          moleculeB.id = B->getVacantID();
        }
    }
  else if(B->getVacantID() == C->getVacantID() ||
          B->getID() == C->getVacantID())
    {
      coordC = coordB;
      if(D)
        {
          if(A->getVacantID() == D->getVacantID() ||
             A->getID() == D->getVacantID())
            {
              coordD = coordA;
            }
          else
            {
              coordD = D->getRandomAdjoiningCoord(coordC, coordC, SearchVacant);
              if(coordD == theNullCoord)
                {
                  return false;
                }
              moleculeA.id = A->getVacantID();
            }
          D->addCoord(coordD);
        }
      else
        {
          //Hard remove the A molecule since it is not used:
          moleculeA.id = A->getVacantID();
        }
    }
  else
    {
      coordC = C->getRandomAdjoiningCoord(coordA, SearchVacant);
      if(coordC == theNullCoord)
        {
          coordC = C->getRandomAdjoiningCoord(coordB, SearchVacant);
          if(coordC == theNullCoord)
            {
              //Only proceed if we can find an adjoining vacant voxel
              //of A or B which can be occupied by C:
              return false;
            }
        }
      if(D)
        {
          coordD = D->getRandomAdjoiningCoord(coordC, coordC, SearchVacant);
          if(coordD == theNullCoord)
            {
              return false;
            }
          D->addCoord(coordD);
        }
      //Hard remove the A molecule since it is not used:
      moleculeA.id = A->getVacantID();
      //Hard remove the B molecule since it is not used:
      moleculeB.id = B->getVacantID();
    }
  C->addCoord(coordC);
  addMoleculeE();
  return true;
}

//zero-coefficient E
//we create a molecule E at random location in the compartment to avoid
//rebinding effect, useful to maintain the concentration of a substrate species
//even after the reaction:
void DiffusionInfluencedReactionProcess::addMoleculeE()
{
  if(!E)
    {
      return;
    } 
  unsigned int coordE(E->getRandomCompCoord(1));
  if(coordE == theNullCoord)
    {
      std::cout << getFullID().asString() << " unable to add molecule E" <<
        std::endl;
    }
  E->addCoord(coordE);
}



void DiffusionInfluencedReactionProcess::finalizeReaction()
{
  //The number of molecules may have changed for both reactant and product
  //species. We need to update SpatiocyteNextReactionProcesses which are
  //dependent on these species:
  for(std::vector<SpatiocyteProcess*>::const_iterator 
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
  if(A->getDimension() == 3 && B->getDimension() == 3)
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
  else if(A->getDimension() != 3 && B->getDimension() != 3)
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
  else if(A->getDimension() == 3 && B->getIsLipid())
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
  else if(A->getIsLipid() && B->getDimension() == 3)
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
  else if(A->getDimension() == 3 && B->getDimension() != 3)
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
  else if(A->getDimension() != 3 && B->getDimension() == 3)
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
                      " [" + getFullID().asString() + 
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


    
      


