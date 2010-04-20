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
  if(theOrder == 1)
    { 
      //nonHD_A -> nonHD_C + nonHD_D:
      if(A && C && D)
        {
          Voxel* moleculeA(A->getRandomMolecule());
          //If the product C is not in the same compartment as A,
          //we need to find a vacant adjoining voxel of A that belongs
          //to the compartment of C:
          Voxel* moleculeC;
          if(A->getCompartment() != C->getCompartment())
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
          //If the product C is not in the same compartment as A,
          //we need to find a vacant adjoining voxel of A that belongs
          //to the compartment of C:
          Voxel* moleculeC;
          if(A->getCompartment() != C->getCompartment())
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
          if(A->getCompartment() != C->getCompartment())
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
        }
      //nonHD_A -> nonHD_C + HD_D:
      else if(A && C && variableD)
        {
        }
      //nonHD_A -> HD_C + nonHD_D:
      else if(A && variableC && D)
        {
        }
      //HD_A -> nonHD_C + nonHD_D:
      else if(variableA && C && D)
        {
        }
      //HD_A -> nonHD_C:
      else if(variableA && C && !D && !variableD)
        {
          Voxel* moleculeC(C->getRandomCompartmentVoxel());
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
      else if(variableA && C && variableD)
        {
        }
      //HD_A -> HD_C + nonHD_D:
      else if(variableA && variableC && D)
        {
        }
    }
  else
    {
      //nonHD_A + HD_B -> nonHD_C + nonHD_D 
      if(A && variableB && C && D)
        {
          Voxel* moleculeA(A->getRandomMolecule());
          //If the product C is not in the same compartment as A,
          //we need to find a vacant adjoining voxel of A that belongs
          //to the compartment of C:
          Voxel* moleculeC;
          if(A->getCompartment() != C->getCompartment())
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
          variableB->addValue(-1);
          A->removeMolecule(moleculeA);
          D->addMolecule(moleculeD);
          C->addMolecule(moleculeC);
        }
      //nonHD_A + HD_B -> nonHD_C
      else if(A && variableB && C && !D)
        {
          Voxel* moleculeA(A->getRandomMolecule());
          //If the product C is not in the same compartment as A,
          //we need to find a vacant adjoining voxel of A that belongs
          //to the compartment of C:
          Voxel* moleculeC;
          if(A->getCompartment() != C->getCompartment())
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
          variableB->addValue(-1);
          A->removeMolecule(moleculeA);
          C->addMolecule(moleculeC);
        }
    }
  ReactionProcess::fire();
}

void SpatiocyteNextReactionProcess::initializeThird()
{
  ReactionProcess::initializeThird();
  if(theOrder == 1) 
    {
      p = k;
    }
  else if(theOrder == 2)
    {
      if(A && variableB)
        {
          Compartment* compB(theSpatiocyteStepper->system2compartment(
                  variableB->getSuperSystem()));
          if(A->getIsVolume() || !compB->isSurface)
            {
              p = k/A->getCompartment()->volume;
              cout << "vol:" << A->getCompartment()->volume << endl; 
            }
          else
            {
              double anArea;
              if(!A->getIsVolume())
                {
                  anArea = A->getCompartment()->area;
                }
              else
                {
                  anArea = compB->area;
                }
              p = k/anArea;
              cout << "area:" << anArea << endl; 
            }
        }
      //A + A -> products
      if(getZeroVariableReferenceOffset() == 1)
        {
          p *= 2;
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

