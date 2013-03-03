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
      reactM = &DiffusionInfluencedReactionProcess::reactWithMultiscaleComp;
    }
  else if(!A->getIsMultiscaleComp() && B->getIsMultiscaleComp())
    {
      N = A;
      M = B;
      isReactWithMultiscaleComp = true;
      reactM = &DiffusionInfluencedReactionProcess::reactWithMultiscaleComp;
    }
  else if(A->getIsMultiscaleComp() && B->getIsMultiscaleComp())
    {
      isReactInMultiscaleComp = true;
      reactM = &DiffusionInfluencedReactionProcess::reactInMultiscaleComp;
    }
  else
    {
      setReactMethod();
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
  unsigned vacantIdx(moleculeM->idx);
  if(M->getIsOnMultiscale())
    {
      vacantIdx = M->getTag(indexM).vacantIdx; 
    }
  if(M_p)
    {
      M_p->addMoleculeInMulti(moleculeM, vacantIdx);
    }
  else
    {
      moleculeM->idx = vacantIdx;
    }
  if(N_p)
    {
      N_p->addMolecule(moleculeN);
    }
  else
    {
      moleculeN->idx = N->getVacantIdx();
    }
}

void DiffusionInfluencedReactionProcess::reactInMultiscaleComp(
                                       Voxel* molA, Voxel* molB,
                                       const unsigned indexA,
                                       const unsigned indexB)
{
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
      if(C->getIsOnMultiscale())
        {
          C->addMoleculeInMulti(molA, vacantIdxA);
        }
      else
        {
          molA->idx = vacantIdxA;
        }
      if(D)
        {
          if(B->isReplaceable(molB, D))
            {
              if(D->getIsOnMultiscale())
                {
                  D->addMoleculeInMulti(molB, vacantIdxB);
                }
              else
                {
                  molB->idx = vacantIdxB;
                }
            }
        }
    }
  else if(B->isReplaceable(molB, C))
    {
      if(C->getIsOnMultiscale())
        {
          C->addMoleculeInMulti(molB, vacantIdxB);
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
                  D->addMoleculeInMulti(molA, vacantIdxA);
                }
              else
                {
                  molA->idx = vacantIdxA;
                }
            }
        }
    }
}

//A + B -> variableC + [D <- molA]
void DiffusionInfluencedReactionProcess::reactVarC_AtoD(
                                                  Voxel* molA, Voxel* molB,
                                                  const unsigned indexA,
                                                  const unsigned indexB)
{
  variableC->addValue(1);
  D->addMolecule(molA, A->getTag(indexA));
  A->softRemoveMolecule(indexA);
  if(A != B)
    {
      B->removeMolecule(indexB);
    }
  else
    {
      //If A == B, indexB is no longer valid after molA is removed,
      //so need to use the current index to remove molB:
      B->removeMolecule(molB);
    }
}

//A + B -> variableC + [D <- molB]
void DiffusionInfluencedReactionProcess::reactVarC_BtoD(
                                                  Voxel* molA, Voxel* molB,
                                                  const unsigned indexA,
                                                  const unsigned indexB)
{
  variableC->addValue(1);
  D->addMolecule(molB, B->getTag(indexB));
  B->softRemoveMolecule(indexB);
  if(A != B)
    {
      A->removeMolecule(indexA);
    }
  else
    {
      A->removeMolecule(molA);
    }
}

//A + B -> variableC + [D <- molN]
void DiffusionInfluencedReactionProcess::reactVarC_NtoD(
                                                  Voxel* molA, Voxel* molB,
                                                  const unsigned indexA,
                                                  const unsigned indexB)
{ 
  Voxel* mol(D->getRandomAdjoiningVoxel(molA, SearchVacant));
  if(!mol)
    {
      mol = D->getRandomAdjoiningVoxel(molB, SearchVacant);
      if(!mol)
        {
          return;
        }
    }
  variableC->addValue(1);
  //TODO: need to use the correct tag here:
  D->addMolecule(mol, A->getTag(indexA));
  A->removeMolecule(indexA);
  if(A != B)
    {
      B->removeMolecule(indexB);
    }
  else
    {
      B->removeMolecule(molB);
    }
}

//A + B -> variableC + [D == molA]
void DiffusionInfluencedReactionProcess::reactVarC_AeqD(
                                                  Voxel* molA, Voxel* molB,
                                                  const unsigned indexA,
                                                  const unsigned indexB)
{
  variableC->addValue(1);
  B->removeMolecule(indexB);
}

//A + B -> variableC + [D == molB]
void DiffusionInfluencedReactionProcess::reactVarC_BeqD(
                                                  Voxel* molA, Voxel* molB,
                                                  const unsigned indexA,
                                                  const unsigned indexB)
{
  variableC->addValue(1);
  A->removeMolecule(indexA);
}

//A + B -> variableD + [C <- molA]
void DiffusionInfluencedReactionProcess::reactVarD_AtoC(
                                                  Voxel* molA, Voxel* molB,
                                                  const unsigned indexA,
                                                  const unsigned indexB)
{
  variableD->addValue(1);
  C->addMolecule(molA, A->getTag(indexA));
  A->softRemoveMolecule(indexA);
  if(A != B)
    {
      B->removeMolecule(indexB);
    }
  else
    {
      B->removeMolecule(molB);
    }
}

//A + B -> variableD + [C <- molB]
void DiffusionInfluencedReactionProcess::reactVarD_BtoC(
                                                  Voxel* molA, Voxel* molB,
                                                  const unsigned indexA,
                                                  const unsigned indexB)
{
  variableD->addValue(1);
  C->addMolecule(molB, B->getTag(indexB));
  B->softRemoveMolecule(indexB);
  if(A != B)
    {
      A->removeMolecule(indexA);
    }
  else
    {
      A->removeMolecule(molA);
    }
}

//A + B -> variableD + [C <- molN]
void DiffusionInfluencedReactionProcess::reactVarD_NtoC(
                                                  Voxel* molA, Voxel* molB,
                                                  const unsigned indexA,
                                                  const unsigned indexB)
{ 
  Voxel* mol(C->getRandomAdjoiningVoxel(molA, SearchVacant));
  if(!mol)
    {
      mol = C->getRandomAdjoiningVoxel(molB, SearchVacant);
      if(!mol)
        {
          return;
        }
    }
  variableD->addValue(1);
  //TODO: need to use the correct tag here:
  C->addMolecule(mol, A->getTag(indexA));
  A->removeMolecule(indexA);
  if(A != B)
    {
      B->removeMolecule(indexB);
    }
  else
    {
      B->removeMolecule(molB);
    }
}


//A + B -> variableD + [C == A]
void DiffusionInfluencedReactionProcess::reactVarD_AeqC(
                                                  Voxel* molA, Voxel* molB,
                                                  const unsigned indexA,
                                                  const unsigned indexB)
{
  variableD->addValue(1);
  B->removeMolecule(indexB);
}

//A + B -> variableD + [C == B]
void DiffusionInfluencedReactionProcess::reactVarD_BeqC(
                                                  Voxel* molA, Voxel* molB,
                                                  const unsigned indexA,
                                                  const unsigned indexB)
{
  variableD->addValue(1);
  A->removeMolecule(indexA);
}


//A + B -> variableC
void DiffusionInfluencedReactionProcess::reactVarC(
                                                  Voxel* molA, Voxel* molB,
                                                  const unsigned indexA,
                                                  const unsigned indexB)
{
  variableC->addValue(1);
  A->removeMolecule(indexA);
  if(A != B)
    {
      B->removeMolecule(indexB);
    }
  else
    {
      B->removeMolecule(molB);
    }
}

//A + B -> [A == C] + [B == D]
void DiffusionInfluencedReactionProcess::reactAeqC_BeqD(
                                                  Voxel* molA, Voxel* molB,
                                                  const unsigned indexA,
                                                  const unsigned indexB)
{
}

//A + B -> [B == C] + [A == D]
void DiffusionInfluencedReactionProcess::reactBeqC_AeqD(
                                                  Voxel* molA, Voxel* molB,
                                                  const unsigned indexA,
                                                  const unsigned indexB)
{
}

//A + B -> [A == C] + [D <- molB]
void DiffusionInfluencedReactionProcess::reactAeqC_BtoD(
                                                  Voxel* molA, Voxel* molB,
                                                  const unsigned indexA,
                                                  const unsigned indexB)
{
  D->addMolecule(molB, B->getTag(indexB));
  B->softRemoveMolecule(indexB);
}

//A + B -> [A == C] + [D <- molN]
void DiffusionInfluencedReactionProcess::reactAeqC_NtoD(
                                                  Voxel* molA, Voxel* molB,
                                                  const unsigned indexA,
                                                  const unsigned indexB)
{
  Voxel* mol(D->getRandomAdjoiningVoxel(molA, SearchVacant));
  if(!mol)
    {
      mol = D->getRandomAdjoiningVoxel(molB, SearchVacant);
      if(!mol)
        {
          return;
        }
    }
  D->addMolecule(mol, B->getTag(indexB));
  B->removeMolecule(indexB);
}

//A + B -> [B == C] + [D <- molA]
void DiffusionInfluencedReactionProcess::reactBeqC_AtoD(
                                                  Voxel* molA, Voxel* molB,
                                                  const unsigned indexA,
                                                  const unsigned indexB)
{
  D->addMolecule(molA, A->getTag(indexA));
  A->softRemoveMolecule(indexA);
}

//A + B -> [B == C] + [D <- molN]
void DiffusionInfluencedReactionProcess::reactBeqC_NtoD(
                                                  Voxel* molA, Voxel* molB,
                                                  const unsigned indexA,
                                                  const unsigned indexB)
{
  Voxel* mol(D->getRandomAdjoiningVoxel(molA, SearchVacant));
  if(!mol)
    {
      mol = D->getRandomAdjoiningVoxel(molB, SearchVacant);
      if(!mol)
        {
          return;
        }
    }
  D->addMolecule(mol, A->getTag(indexA));
  A->removeMolecule(indexA);
}

//A + B -> [C <- molB] + [A == D]
void DiffusionInfluencedReactionProcess::reactBtoC_AeqD(
                                                  Voxel* molA, Voxel* molB,
                                                  const unsigned indexA,
                                                  const unsigned indexB)
{
  C->addMolecule(molB, B->getTag(indexB));
  B->softRemoveMolecule(indexB);
}

//A + B -> [C <- molN] + [A == D]
void DiffusionInfluencedReactionProcess::reactNtoC_AeqD(
                                                  Voxel* molA, Voxel* molB,
                                                  const unsigned indexA,
                                                  const unsigned indexB)
{
  Voxel* mol(C->getRandomAdjoiningVoxel(molA, SearchVacant));
  if(!mol)
    {
      mol = C->getRandomAdjoiningVoxel(molB, SearchVacant);
      if(!mol)
        {
          return;
        }
    }
  C->addMolecule(mol, B->getTag(indexB));
  B->removeMolecule(indexB);
}

//A + B -> [C <- molA] + [B == D]
void DiffusionInfluencedReactionProcess::reactAtoC_BeqD(
                                                  Voxel* molA, Voxel* molB,
                                                  const unsigned indexA,
                                                  const unsigned indexB)
{
  C->addMolecule(molA, A->getTag(indexA));
  A->softRemoveMolecule(indexA);
}

//A + B -> [C <- molN] + [B == D]
void DiffusionInfluencedReactionProcess::reactNtoC_BeqD(
                                                  Voxel* molA, Voxel* molB,
                                                  const unsigned indexA,
                                                  const unsigned indexB)
{
  Voxel* mol(C->getRandomAdjoiningVoxel(molA, SearchVacant));
  if(!mol)
    {
      mol = C->getRandomAdjoiningVoxel(molB, SearchVacant);
      if(!mol)
        {
          return;
        }
    }
  C->addMolecule(mol, A->getTag(indexA));
  A->removeMolecule(indexA);
}

//A + B -> [C <- molA] + [D <- molB]
void DiffusionInfluencedReactionProcess::reactAtoC_BtoD(
                                                  Voxel* molA, Voxel* molB,
                                                  const unsigned indexA,
                                                  const unsigned indexB)
{
  C->addMolecule(molA, A->getTag(indexA));
  A->softRemoveMolecule(indexA);
  if(A != B)
    {
      D->addMolecule(molB, B->getTag(indexB));
      B->softRemoveMolecule(indexB);
    }
  else
    {
      Tag& aTag(B->getTag(molB));
      B->softRemoveMolecule(molB);
      D->addMolecule(molB, aTag);
    }
}

//A + B -> [C <- molA] + [D <- molN]
void DiffusionInfluencedReactionProcess::reactAtoC_NtoD(
                                                  Voxel* molA, Voxel* molB,
                                                  const unsigned indexA,
                                                  const unsigned indexB)
{
  Voxel* mol(D->getRandomAdjoiningVoxel(molA, SearchVacant));
  if(!mol)
    {
      mol = D->getRandomAdjoiningVoxel(molB, SearchVacant);
      if(!mol)
        {
          return;
        }
    }
  D->addMolecule(mol, B->getTag(indexB));
  C->addMolecule(molA, A->getTag(indexA));
  A->softRemoveMolecule(indexA);
  if(A != B)
    {
      B->removeMolecule(indexB);
    }
  else
    {
      B->removeMolecule(molB);
    }
}

//A + B -> [C <- molB] + [D <- molA]
void DiffusionInfluencedReactionProcess::reactBtoC_AtoD(
                                                  Voxel* molA, Voxel* molB,
                                                  const unsigned indexA,
                                                  const unsigned indexB)
{
  C->addMolecule(molB, B->getTag(indexB));
  B->softRemoveMolecule(indexB);
  if(A != B)
    { 
      D->addMolecule(molA, A->getTag(indexA));
      A->softRemoveMolecule(indexA);
    }
  else
    {
      Tag& aTag(A->getTag(molA));
      A->softRemoveMolecule(molA);
      D->addMolecule(molA, aTag);
    }
}

//A + B -> [C <- molB] + [D <- molN]
void DiffusionInfluencedReactionProcess::reactBtoC_NtoD(
                                                  Voxel* molA, Voxel* molB,
                                                  const unsigned indexA,
                                                  const unsigned indexB)
{
  Voxel* mol(D->getRandomAdjoiningVoxel(molA, SearchVacant));
  if(!mol)
    {
      mol = D->getRandomAdjoiningVoxel(molB, SearchVacant);
      if(!mol)
        {
          return;
        }
    }
  D->addMolecule(mol, A->getTag(indexA));
  C->addMolecule(molB, B->getTag(indexB));
  B->softRemoveMolecule(indexB);
  if(A != B)
    {
      A->removeMolecule(indexA);
    }
  else
    {
      A->removeMolecule(molA);
    }
}

//A + B -> [C <- molN] + [D <- molN]
void DiffusionInfluencedReactionProcess::reactNtoC_NtoD(
                                                  Voxel* molA, Voxel* molB,
                                                  const unsigned indexA,
                                                  const unsigned indexB)
{
  Voxel* molC(C->getRandomAdjoiningVoxel(molA, SearchVacant));
  if(!molC)
    {
      molC = C->getRandomAdjoiningVoxel(molB, SearchVacant);
      if(!molC)
        {
          return;
        }
    }
  Voxel* molD(D->getRandomAdjoiningVoxel(molA, molC, SearchVacant));
  if(!molD)
    {
      molD = D->getRandomAdjoiningVoxel(molB, molC, SearchVacant);
      if(!molD)
        {
          return;
        }
    }
  C->addMolecule(molC, A->getTag(indexA));
  D->addMolecule(molD, B->getTag(indexB));
  A->removeMolecule(indexA);
  if(A != B)
    {
      B->removeMolecule(indexB);
    }
  else
    {
      B->removeMolecule(molB);
    }
}


//A + B -> [A == C]
void DiffusionInfluencedReactionProcess::reactAeqC(
                                                  Voxel* molA, Voxel* molB,
                                                  const unsigned indexA,
                                                  const unsigned indexB)
{
  B->removeMolecule(indexB);
}

//A + B -> [B == C]
void DiffusionInfluencedReactionProcess::reactBeqC(
                                                  Voxel* molA, Voxel* molB,
                                                  const unsigned indexA,
                                                  const unsigned indexB)
{
  A->removeMolecule(indexA);
}

//A + B -> [C <- molA]
void DiffusionInfluencedReactionProcess::reactAtoC(
                                                  Voxel* molA, Voxel* molB,
                                                  const unsigned indexA,
                                                  const unsigned indexB)
{
  C->addMolecule(molA, A->getTag(indexA));
  A->softRemoveMolecule(indexA);
  if(A != B)
    {
      B->removeMolecule(indexB);
    }
  else
    {
      B->removeMolecule(molB);
    }
}

//A + B -> [C <- molB]
void DiffusionInfluencedReactionProcess::reactBtoC(
                                                  Voxel* molA, Voxel* molB,
                                                  const unsigned indexA,
                                                  const unsigned indexB)
{
  C->addMolecule(molB, B->getTag(indexB));
  B->softRemoveMolecule(indexB);
  if(A != B)
    {
      A->removeMolecule(indexA);
    }
  else
    {
      A->removeMolecule(molA);
    }
}

//A + B -> [C <- molN]
void DiffusionInfluencedReactionProcess::reactNtoC(
                                                  Voxel* molA, Voxel* molB,
                                                  const unsigned indexA,
                                                  const unsigned indexB)
{
  Voxel* mol(C->getRandomAdjoiningVoxel(molA, SearchVacant));
  if(!mol)
    {
      mol = C->getRandomAdjoiningVoxel(molB, SearchVacant);
      if(!mol)
        {
          return;
        }
    }
  C->addMolecule(mol, A->getTag(indexA));
  A->removeMolecule(indexA);
  if(A != B)
    {
      B->removeMolecule(indexB);
    }
  else
    {
      B->removeMolecule(molB);
    }
}

void DiffusionInfluencedReactionProcess::setReactMethod()
{
  if(variableC && D)
    {
      if(A == D)
        {
          //A + B -> variableC + [A == D]
          reactM = &DiffusionInfluencedReactionProcess::reactVarC_AeqD;
        }
      else if(B == D)
        {
          //A + B -> variableC + [B == D]
          reactM = &DiffusionInfluencedReactionProcess::reactVarC_BeqD;
        }
      else
        { 
          if(A->isReplaceable(D))
            {
              //A + B -> variableC + [D <- molA]
              reactM = &DiffusionInfluencedReactionProcess::reactVarC_AtoD;
            }
          else if(B->isReplaceable(D))
            {
              //A + B -> variableC + [D <- molB]
              reactM = &DiffusionInfluencedReactionProcess::reactVarC_BtoD;
            }
          else
            {
              //A + B -> variableC + [D <- molN]
              reactM = &DiffusionInfluencedReactionProcess::reactVarC_NtoD;
            }
        }
    }
  else if(variableD && C)
    {
      if(A == C)
        {
          //A + B -> variableD + [A == C]
          reactM = &DiffusionInfluencedReactionProcess::reactVarD_AeqC;
        }
      else if(B == C)
        {
          //A + B -> variableD + [B == C]
          reactM = &DiffusionInfluencedReactionProcess::reactVarD_BeqC;
        }
      else
        { 
          if(A->isReplaceable(C))
            {
              //A + B -> variableD + [C <- molA]
              reactM = &DiffusionInfluencedReactionProcess::reactVarD_AtoC;
            }
          else if(B->isReplaceable(C))
            {
              //A + B -> variableD + [C <- molB]
              reactM = &DiffusionInfluencedReactionProcess::reactVarD_BtoC;
            }
          else
            {
              //A + B -> variableD + [C <- molN]
              reactM = &DiffusionInfluencedReactionProcess::reactVarD_NtoC;
            }
        }
    }
  else if(variableC)
    {
      //A + B -> variableC
      reactM = &DiffusionInfluencedReactionProcess::reactVarC;
    }
  else if(D)
    {
      if(A == C && B == D)
        {
          //A + B -> [A == C] + [B == D]
          reactM = &DiffusionInfluencedReactionProcess::reactAeqC_BeqD;
        }
      else if(B == C && A == D)
        {
          //A + B -> [B == C] + [A == D]
          reactM = &DiffusionInfluencedReactionProcess::reactBeqC_AeqD;
        }
      else if(A == C)
        {
          if(B->isReplaceable(D))
            {
              //A + B -> [A == C] + [D <- molB]
              reactM = &DiffusionInfluencedReactionProcess::reactAeqC_BtoD;
            }
          else
            {
              //A + B -> [A == C] + [D <- molN]
              reactM = &DiffusionInfluencedReactionProcess::reactAeqC_NtoD;
            }
        }
      else if(B == C)
        {
          if(A->isReplaceable(D))
            {
              //A + B -> [B == C] + [D <- molA]
              reactM = &DiffusionInfluencedReactionProcess::reactBeqC_AtoD;
            }
          else
            {
              //A + B -> [B == C] + [D <- molN]
              reactM = &DiffusionInfluencedReactionProcess::reactBeqC_NtoD;
            }
        }
      else if(A == D)
        {
          if(B->isReplaceable(C))
            {
              //A + B -> [C <- molB] + [A == D]
              reactM = &DiffusionInfluencedReactionProcess::reactBtoC_AeqD;
            }
          else
            {
              //A + B -> [C <- molN] + [A == D]
              reactM = &DiffusionInfluencedReactionProcess::reactNtoC_AeqD;
            }
        }
      else if(B == D)
        {
          if(A->isReplaceable(C))
            {
              //A + B -> [C <- molA] + [B == D]
              reactM = &DiffusionInfluencedReactionProcess::reactAtoC_BeqD;
            }
          else
            {
              //A + B -> [C <- molN] + [B == D]
              reactM = &DiffusionInfluencedReactionProcess::reactNtoC_BeqD;
            }
        }
      else
        {
          if(A->isReplaceable(C))
            {
              if(B->isReplaceable(D))
                {
                  //A + B -> [C <- molA] + [D <- molB]
                  reactM = &DiffusionInfluencedReactionProcess::reactAtoC_BtoD;
                }
              else
                {
                  //A + B -> [C <- molA] + [D <- molN]
                  reactM = &DiffusionInfluencedReactionProcess::reactAtoC_NtoD;
                }
            }
          else if(B->isReplaceable(C))
            {
              if(A->isReplaceable(D))
                {
                  //A + B -> [C <- molB] + [D <- molA]
                  reactM = &DiffusionInfluencedReactionProcess::reactBtoC_AtoD;
                }
              else
                {
                  //A + B -> [C <- molB] + [D <- molN]
                  reactM = &DiffusionInfluencedReactionProcess::reactBtoC_NtoD;
                }
            }
          else
            {
              //A + B -> [C <- molN] + [D <- molN]
              reactM = &DiffusionInfluencedReactionProcess::reactNtoC_NtoD;
            }
        }
    }
  else
    {
      if(A == C)
        {
          //A + B -> [A == C]
          reactM = &DiffusionInfluencedReactionProcess::reactAeqC;
        }
      else if(B == C)
        {
          //A + B -> [B == C]
          reactM = &DiffusionInfluencedReactionProcess::reactBeqC;
        }
      else if(A->isReplaceable(C))
        {
          //A + B -> [C <- molA]
          reactM = &DiffusionInfluencedReactionProcess::reactAtoC;
        }
      else if(B->isReplaceable(C))
        {
          //A + B -> [C <- molB]
          reactM = &DiffusionInfluencedReactionProcess::reactBtoC;
        }
      else
        {
          //A + B -> [C <- molN]
          reactM = &DiffusionInfluencedReactionProcess::reactNtoC;
        }
    }
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


    
      


