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

#include "SpatiocyteNextReactionProcess.hpp"
#include "SpatiocyteSpecies.hpp"

LIBECS_DM_INIT(SpatiocyteNextReactionProcess, Process);

void SpatiocyteNextReactionProcess::fire()
{
  if(theOrder == 0)
    {
      if(C)
        { 
          Voxel* moleculeC(C->getRandomCompVoxel());
          if(moleculeC == NULL)
            {
              requeue();
              return;
            }
          C->addMolecule(moleculeC);
        }
      else if(variableC)
        {
          variableC->addValue(1);
        }
    }
  else if(theOrder == 1)
    { 
      //nonHD_A -> nonHD_C + nonHD_D:
      if(A && C && D)
        {
          Voxel* moleculeA(A->getRandomMolecule());
          //If the product C is not in the same Comp as A,
          //we need to find a vacant adjoining voxel of A that belongs
          //to the Comp of C:
          Voxel* moleculeC;
          if(A->getComp() != C->getComp())
            {
              moleculeC = C->getRandomAdjoiningVoxel(moleculeA);
              //Only proceed if we can find an adjoining vacant voxel
              //of A which can be occupied by C:
              if(moleculeC == NULL)
                {
                  requeue();
                  return;
                }
            }
          else
            {
              moleculeC = moleculeA;
            }
          Voxel* moleculeD(D->getRandomAdjoiningVoxel(moleculeC, moleculeC));
          //Only proceed if we can find an adjoining vacant voxel
          //of A which can be occupied by D:
          if(moleculeD == NULL)
            {
              requeue();
              return;
            }
          D->addMolecule(moleculeD);
          A->removeMolecule(moleculeA);
          C->addMolecule(moleculeC);
        }
      //nonHD_A -> nonHD_C:
      else if(A && C && !D && !variableD)
        {
          Voxel* moleculeA(A->getRandomMolecule());
          //If the product C is not in the same Comp as A,
          //we need to find a vacant adjoining voxel of A that belongs
          //to the Comp of C:
          Voxel* moleculeC;
          if(A->getComp() != C->getComp())
            {
              moleculeC = C->getRandomAdjoiningVoxel(moleculeA);
              //Only proceed if we can find an adjoining vacant voxel
              //of A which can be occupied by C:
              if(moleculeC == NULL)
                {
                  requeue();
                  return;
                }
            }
          else
            {
              moleculeC = moleculeA;
            }
          A->removeMolecule(moleculeA);
          C->addMolecule(moleculeC);
        }
      //nonHD_A -> HD_C + HD_D:
      else if(A && variableC && variableD && !variableE)
        {
          Voxel* moleculeA(A->getRandomMolecule());
          A->removeMolecule(moleculeA);
          variableC->addValue(1);
          variableD->addValue(1);
        }
      //nonHD_A -> nonHD_C + HD_D + HD_E:
      else if(A && C && variableD && variableE)
        {
          Voxel* moleculeA(A->getRandomMolecule());
          Voxel* moleculeC;
          if(A->getComp() != C->getComp())
            {
              moleculeC = C->getRandomAdjoiningVoxel(moleculeA);
              //Only proceed if we can find an adjoining vacant voxel
              //of A which can be occupied by C:
              if(moleculeC == NULL)
                {
                  requeue();
                  return;
                }
            }
          else
            {
              moleculeC = moleculeA;
            }
          A->removeMolecule(moleculeA);
          C->addMolecule(moleculeC);
          variableD->addValue(1);
          variableE->addValue(1);
        }
      //nonHD_A -> HD_C:
      else if(A && variableC && !D && !variableD)
        {
          Voxel* moleculeA(A->getRandomMolecule());
          A->removeMolecule(moleculeA);
          variableC->addValue(1);
        }
      //nonHD_A -> nonHD_C + HD_D:
      //nonHD_A -> HD_C + nonHD_D:
      else if(A && ((variableC && D) || (C && variableD)))
        {
          Variable* HD_p(variableC);
          Species* nonHD_p(D);
          if(variableD)
            {
              HD_p = variableD;
              nonHD_p = C;
            }
          Voxel* moleculeA(A->getRandomMolecule());
          Voxel* molecule;
          if(A->getComp() != nonHD_p->getComp())
            {
              molecule = nonHD_p->getRandomAdjoiningVoxel(moleculeA);
              //Only proceed if we can find an adjoining vacant voxel
              //of A which can be occupied by nonHD:
              if(molecule == NULL)
                {
                  requeue();
                  return;
                }
            }
          else
            {
              molecule = moleculeA;
            }
          A->removeMolecule(moleculeA);
          nonHD_p->addMolecule(molecule);
          HD_p->addValue(1);
        }
      //HD_A -> nonHD_C + nonHD_D:
      else if(variableA && C && D)
        {
        }
      //HD_A -> nonHD_C:
      else if(variableA && C && !D && !variableD)
        {
          Voxel* moleculeC(C->getRandomCompVoxel());
          if(moleculeC == NULL)
            {
              requeue();
              return;
            }
          variableA->addValue(-1);
          C->addMolecule(moleculeC);
        }
      //HD_A -> HD_C + HD_D:
      else if(variableA && variableC && variableD)
        {
          variableA->addValue(-1);
          variableC->addValue(1);
          variableD->addValue(1);
        }
      //HD_A -> HD_C:
      else if(variableA && variableC && !D && !variableD)
        {
          variableA->addValue(-1);
          variableC->addValue(1);
        }
      //HD_A -> nonHD_C + HD_D:
      //HD_A -> HD_C + nonHD_D:
      else if(variableA && ((variableC && D) || (C && variableD)))
        {
          Variable* HD_p(variableC);
          Species* nonHD_p(D);
          if(variableD)
            {
              HD_p = variableD;
              nonHD_p = C;
            }
          Voxel* molecule(nonHD_p->getRandomCompVoxel());
          if(molecule == NULL)
            {
              requeue();
              return;
            }
          variableA->addValue(-1);
          nonHD_p->addMolecule(molecule);
          HD_p->addValue(1);
        }
    }
  else
    {
      Species* nonHD(A);
      Variable* HD(variableB);
      if(B)
        {
          nonHD = B;
          HD = variableA;
        }
      //nonHD + HD -> nonHD + nonHD: 
      //HD + nonHD -> nonHD + nonHD: 
      if(C && D)
        {
          Voxel* moleculeNonHD(nonHD->getRandomMolecule());
          //If the product C is not in the same Comp as nonHD,
          //we need to find a vacant adjoining voxel of nonHD that belongs
          //to the Comp of C:
          Voxel* moleculeC;
          if(nonHD->getComp() != C->getComp())
            {
              moleculeC = C->getRandomAdjoiningVoxel(moleculeNonHD);
              //Only proceed if we can find an adjoining vacant voxel
              //of nonHD which can be occupied by C:
              if(moleculeC == NULL)
                {
                  requeue();
                  return;
                }
            }
          else
            {
              moleculeC = moleculeNonHD;
            }
          Voxel* moleculeD(D->getRandomAdjoiningVoxel(moleculeC, moleculeC));
          //Only proceed if we can find an adjoining vacant voxel
          //of A which can be occupied by D:
          if(moleculeD == NULL)
            {
              requeue();
              return;
            }
          HD->addValue(-1);
          nonHD->removeMolecule(moleculeNonHD);
          D->addMolecule(moleculeD);
          C->addMolecule(moleculeC);
        }
      //nonHD + HD -> nonHD:
      //HD + nonHD -> nonHD:
      else if(C && !D)
        {
          Voxel* moleculeNonHD(nonHD->getRandomMolecule());
          //If the product C is not in the same Comp as nonHD,
          //we need to find a vacant adjoining voxel of nonHD that belongs
          //to the Comp of C:
          Voxel* moleculeC;
          if(nonHD->getComp() != C->getComp())
            {
              moleculeC = C->getRandomAdjoiningVoxel(moleculeNonHD);
              //Only proceed if we can find an adjoining vacant voxel
              //of nonHD which can be occupied by C:
              if(moleculeC == NULL)
                {
                  requeue();
                  return;
                }
            }
          else
            {
              moleculeC = moleculeNonHD;
            }
          HD->addValue(-1);
          nonHD->removeMolecule(moleculeNonHD);
          C->addMolecule(moleculeC);
        }
    }
  ReactionProcess::fire();
}

void SpatiocyteNextReactionProcess::initializeThird()
{
  ReactionProcess::initializeThird();
  if(p != -1)
    {
      return;
    }
  if(theOrder == 0)
    {
      Comp* compC(NULL);
      if(C)
        {
          compC = C->getComp();
        }
      else
        {
          compC = theSpatiocyteStepper->system2Comp(
                             variableC->getSuperSystem());
        }
      double aSpace(0);
      if(compC->isSurface)
        {
          aSpace = compC->actualArea;
        }
      else
        {
          aSpace = compC->actualVolume;
        }
      p = k*aSpace;
    }
  else if(theOrder == 1) 
    {
      p = k;
    }
  else if(theOrder == 2)
    {
      Comp* compA(NULL);
      Comp* compB(NULL);
      Comp* compC(NULL);
      Comp* compD(NULL);
      if(A)
        {
          compA = A->getComp();
        }
      else
        {
          compA = theSpatiocyteStepper->system2Comp(
                             variableA->getSuperSystem());
        }
      if(B)
        {
          compB = B->getComp();
        }
      else
        {
          compB = theSpatiocyteStepper->system2Comp(
                             variableB->getSuperSystem());
        }
      if(C)
        {
          compC = C->getComp();
        }
      else
        {
          compC = theSpatiocyteStepper->system2Comp(
                             variableC->getSuperSystem());
        }
      if(D)
        {
          compD = D->getComp();
        }
      else if(variableD)
        {
          compD = theSpatiocyteStepper->system2Comp(
                             variableD->getSuperSystem());
        }
      //If there are two products, both C and D must belong to the
      //same Comp:
      if(compD && compD != compC)
        {
          NEVER_GET_HERE;
        }
      double aVolume(compA->actualVolume);
      double anArea(compA->actualArea);
      if(compB->isSurface)
        {
          anArea = compB->actualArea;
        }
      else
        {
          aVolume = compB->actualVolume;
        }
      if(compC->isSurface)
        {
          anArea = compC->actualArea;
        }
      else
        {
          aVolume = compC->actualVolume;
        }
      //If volume (+volume) = k(volume)(volume) or
      //   surface (+surface) = k(volume)(surface) or
      //   surface (+surface) = k(surface)(volume)
      if((!compC->isSurface && !compA->isSurface && !compB->isSurface) ||
         (compC->isSurface && !compA->isSurface && compB->isSurface) ||
         (compC->isSurface && compA->isSurface && !compB->isSurface))
        {
          p = k/aVolume;
        }
      //If surface (+surface) = k(surface)(surface) or
      //   volume (+volume) = k(volume)(surface) or
      //   volume (+volume) = k(surface)(volume)
      else if((compC->isSurface && compA->isSurface && compB->isSurface) ||
              (!compC->isSurface && !compA->isSurface && compB->isSurface) ||
              (!compC->isSurface && compA->isSurface && !compB->isSurface))
        {
          p = k/anArea;
        }
      else
        {
          NEVER_GET_HERE;
        }
      //A + A -> products
      if(getZeroVariableReferenceOffset() == 1)
        {
          p = k;
        }
    }
  else
    {
      NEVER_GET_HERE;
    } 
}


GET_METHOD_DEF(Real, StepInterval, SpatiocyteNextReactionProcess)
{
  return getPropensity_R()*
    (-log(gsl_rng_uniform_pos(getStepper()->getRng())));
}

