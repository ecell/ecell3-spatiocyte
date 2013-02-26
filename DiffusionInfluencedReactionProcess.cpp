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
		    "]: This process cannot have a HD substrate species: " + 
        getIDString(variableA));
    }
  if(variableB)
    {
      THROW_EXCEPTION(ValueError, String(
        getPropertyInterface().getClassName()) + " [" + getFullID().asString() +
		    "]: This process cannot have a HD substrate species: " + 
        getIDString(variableB));
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
  if(A->getIsMultiscaleComp() && !B->getIsMultiscaleComp())
    {
      M = A;
      N = B;
      isReactWithMultiscaleComp = true;
    }
  else if(!A->getIsMultiscaleComp() && B->getIsMultiscaleComp())
    {
      N = A;
      M = B;
      isReactWithMultiscaleComp = true;
    }
  else if(A->getIsMultiscaleComp() && B->getIsMultiscaleComp())
    {
      isReactInMultiscaleComp = true;
    }
  if(isReactWithMultiscaleComp)
    {
      if(C->getIsMultiscaleComp())
        {
          M_p = C;
          if(D)
            {
              N_p = D;
            }
        }
      else
        {
          N_p = C;
          if(D)
            {
              M_p = D;
            }
        }
    }
}


//M -> isMultiscaleComp (isMultiscale or isOnMultiscale)
//N -> normal species (not isMultiscaleComp)
void DiffusionInfluencedReactionProcess::reactWithMultiscaleComp(
                                       Voxel* moleculeN, Voxel* moleculeM,
                                       const unsigned indexN,
                                       const unsigned indexM)
{
  std::cout << "checking before" << std::endl;
  theSpatiocyteStepper->checkSpecies();
  unsigned vacantIdx(moleculeM->idx);
  std::cout << "vacIdx:" << vacantIdx << std::endl;
  if(M->getIsOnMultiscale())
    {
      vacantIdx = M->getTag(indexM).vacantIdx; 
      std::cout << "b vacIdx:" << vacantIdx << " sp:" << M->getIDString(vacantIdx/theStride) << std::endl;
    }
  if(M_p)
    {
      std::cout << "1 M:" << M->getIDString() << " N:" << N->getIDString() << 
        " M_p:" << M_p->getIDString() << " vacIdx:" << vacantIdx << " sp:" << M->getIDString(vacantIdx/theStride) << std::endl;
      M_p->addMolecule(moleculeM, vacantIdx);
      for(unsigned i(0); i != M_p->size(); ++i)
        {
          if(!M_p->getTag(i).vacantIdx)
            {
              std::cout << "error in vacantIdx:" << i << std::endl;
            }
        }
      std::cout << "1 M added sp:" << M_p->getIDString(moleculeM->idx/theStride) << std::endl;
    }
  else
    {
      std::cout << "2 M:" << M->getIDString() << " N:" << N->getIDString() << std::endl;
      moleculeM->idx = vacantIdx;
    }
  if(N_p)
    {
      std::cout << "3 M:" << M->getIDString() << " N:" << N->getIDString() <<
        " N_p:" <<  N_p->getIDString() << std::endl;
      N_p->addMolecule(moleculeN);
    }
  else
    {
      std::cout << "4 M:" << M->getIDString() << " N:" << N->getIDString() << " vacIdx:" << vacantIdx << " vac:" << N->getIDString(N->getVacantIdx()/theStride) << std::endl;
      moleculeN->idx = N->getVacantIdx();
    }
}

void DiffusionInfluencedReactionProcess::reactInMultiscaleComp(
                                       Voxel* molA, Voxel* molB,
                                       const unsigned indexA,
                                       const unsigned indexB)
{
  std::cout << "in multi checking before" << std::endl;
  theSpatiocyteStepper->checkSpecies();
  unsigned vacantIdxA(molA->idx);
  if(A->getIsOnMultiscale())
    {
      vacantIdxA = A->getTag(indexA).vacantIdx; 
    }
  unsigned vacantIdxB(molB->idx);
  if(B->getIsOnMultiscale())
    {
      vacantIdxB = B->getTag(indexB).vacantIdx; 
    }
  if(A->isReplaceable(molA, C))
    {
      std::cout << "1 A:" << A->getIDString() << " B:" << B->getIDString() << " vacIdxA:" << vacantIdxA << " vacB:" << vacantIdxB << " C:" << C->getIDString() << std::endl;
      if(C->getIsOnMultiscale())
        {
          std::cout << "a" << std::endl;
          C->addMolecule(molA, vacantIdxA);
        }
      else
        {
          std::cout << "b" << std::endl;
          molA->idx = vacantIdxA;
        }
      if(D)
        {
          if(B->isReplaceable(molB, D))
            {
              if(D->getIsOnMultiscale())
                {
          std::cout << "c" << std::endl;
                  D->addMolecule(molB, vacantIdxB);
                }
              else
                {
          std::cout << "d" << std::endl;
                  molB->idx = vacantIdxB;
                }
            }
        }
    }
  else if(B->isReplaceable(molB, C))
    {
      std::cout << "2 A:" << A->getIDString() << " B:" << B->getIDString() << " vacIdxA:" << vacantIdxA << " vacB:" << vacantIdxB << " C:" << C->getIDString() << std::endl;
      if(C->getIsOnMultiscale())
        {
          C->addMolecule(molB, vacantIdxB);
        }
      else
        {
          molB->idx = vacantIdxB;
        }
      if(D)
        {
          if(A->isReplaceable(molA, D))
            {
              if(D->getIsOnMultiscale())
                {
                  D->addMolecule(molA, vacantIdxA);
                }
              else
                {
                  molA->idx = vacantIdxA;
                }
            }
        }
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
bool DiffusionInfluencedReactionProcess::react(Voxel* molA, Voxel* molB,
                                               const unsigned indexA,
                                               const unsigned indexB)
{
  moleculeA = molA;
  moleculeB = molB;
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
      if(A->isReplaceable(moleculeA, nonHD_p))
        {
          moleculeP = moleculeA;
          //Hard remove the B molecule, since nonHD_p is in a different Comp:
          moleculeB->idx = B->getVacantIdx();
        }
      else if(B->isReplaceable(moleculeB, nonHD_p))
        {
          moleculeP = moleculeB;
          //Hard remove the A molecule, since nonHD_p is in a different Comp:
          moleculeA->idx = A->getVacantIdx();
        }
      else
        { 
          moleculeP = nonHD_p->getRandomAdjoiningVoxel(moleculeA, SearchVacant);
          //Only proceed if we can find an adjoining vacant voxel
          //of A which can be occupied by C:
          if(moleculeP == NULL)
            {
              moleculeP = nonHD_p->getRandomAdjoiningVoxel(moleculeB,
                                                           SearchVacant);
              if(moleculeP == NULL)
                {
                  return false;
                }
            }
          //Hard remove the A molecule, since nonHD_p is in a different Comp:
          moleculeA->idx = A->getVacantIdx();
          //Hard remove the B molecule, since nonHD_p is in a different Comp:
          moleculeB->idx = B->getVacantIdx();
        }
      HD_p->addValue(1);
      nonHD_p->addMolecule(moleculeP, A->getTag(indexA));
      return true;
    }
  //nonHD_A + nonHD_B -> HD_C:
  else if(variableC && !D && !variableD)
    {

      //Hard remove the A molecule, since nonHD_p is in a different Comp:
      moleculeA->idx = A->getVacantIdx();
      //Hard remove the B molecule, since nonHD_p is in a different Comp:
      moleculeB->idx = B->getVacantIdx();
      variableC->addValue(1);
      return true;
    }

  if(A->isReplaceable(moleculeA, C))
    {
      moleculeC = moleculeA;
      if(D)
        {
          if(B->isReplaceable(moleculeB, D))
            {
              moleculeD = moleculeB;
            }
          else
            {
              moleculeD = D->getRandomAdjoiningVoxel(moleculeC, moleculeC,
                                                     SearchVacant);
              if(moleculeD == NULL)
                {
                  return false;
                }
              moleculeB->idx = B->getVacantIdx();
            }
          D->addMolecule(moleculeD, B->getTag(indexB));
        }
      else
        {
          //Hard remove the B molecule since it is not used:
          moleculeB->idx = B->getVacantIdx();
        }
    }
  else if(B->isReplaceable(moleculeB, C))
    {
      moleculeC = moleculeB;
      if(D)
        {
          if(A->isReplaceable(moleculeA, D))
            {
              moleculeD = moleculeA;
            }
          else
            {
              moleculeD = D->getRandomAdjoiningVoxel(moleculeC, moleculeC,
                                                     SearchVacant);
              if(moleculeD == NULL)
                {
                  return false;
                }
              moleculeA->idx = A->getVacantIdx();
            }
          D->addMolecule(moleculeD, B->getTag(indexB));
        }
      else
        {
          //Hard remove the A molecule since it is not used:
          moleculeA->idx = A->getVacantIdx();
        }
    }
  else
    {
      moleculeC = C->getRandomAdjoiningVoxel(moleculeA, SearchVacant);
      if(moleculeC == NULL)
        {
          moleculeC = C->getRandomAdjoiningVoxel(moleculeB, SearchVacant);
          if(moleculeC == NULL)
            {
              //Only proceed if we can find an adjoining vacant voxel
              //of A or B which can be occupied by C:
              return false;
            }
        }
      if(D)
        {
          moleculeD = D->getRandomAdjoiningVoxel(moleculeC, moleculeC,
                                                 SearchVacant);
          if(moleculeD == NULL)
            {
              return false;
            }
          D->addMolecule(moleculeD, B->getTag(indexB));
        }
      //Hard remove the A molecule since it is not used:
      moleculeA->idx = A->getVacantIdx();
      //Hard remove the B molecule since it is not used:
      moleculeB->idx = B->getVacantIdx();
    }
  C->addMolecule(moleculeC, A->getTag(indexA));
  addMoleculeE();
  addMoleculeF();
  return true;
}

//positive-coefficient F
void DiffusionInfluencedReactionProcess::addMoleculeF()
{
  if(!F)
    {
      return;
    }
  moleculeF = F->getRandomAdjoiningVoxel(moleculeC, SearchVacant);
  if(moleculeF == NULL)
    {
      moleculeF = F->getRandomAdjoiningVoxel(moleculeD, SearchVacant);
      if(moleculeF == NULL)
        {
          return;
        }
    }
  F->addMolecule(moleculeF);
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
  moleculeE = E->getRandomCompVoxel(1);
  if(moleculeE == NULL)
    {
      std::cout << getFullID().asString() << " unable to add molecule E" <<
        std::endl;
      return;
    }
  E->addMolecule(moleculeE);
}

void DiffusionInfluencedReactionProcess::finalizeReaction()
{
  //The number of molecules may have changed for both reactant and product
  //species. We need to update SpatiocyteNextReactionProcesses which are
  //dependent on these species:
  for(std::vector<SpatiocyteProcess*>::const_iterator 
      i(theInterruptedProcesses.begin());
      i!=theInterruptedProcesses.end(); ++i)
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


    
      


