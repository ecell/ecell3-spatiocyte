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
  /*
  std::cout << std::endl << "before:" << getFullID().asString() << std::endl;
  Real aValue1(theVariableReferenceVector[0].getVariable()->getValue());
  Real aValue2(theVariableReferenceVector[1].getVariable()->getValue());
  std::cout << "aValue1:" << aValue1 << " " << variableA->getValue();
  std::cout << " aValue2:" << aValue2 << " " << B->size() << std::endl;
  */

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
          if(!reactACD(A, C, D))
            {
              return;
            }
        }
      //nonHD_A -> nonHD_C:
      else if(A && C && !D && !variableD)
        {
          if(BindingSite == -1)
            {
              if(!reactAC(A, C))
                {
                  return;
                }
            }
          else if(!reactACbind(A, C))
            {
              return;
            }
        }
      //nonHD_A -> HD_C + HD_D:
      else if(A && variableC && variableD)
        {
          Voxel* moleculeA(A->getRandomMolecule());
          A->removeMolecule(moleculeA);
          variableC->addValue(1);
          variableD->addValue(1);
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
          if(reactAC(A, nonHD_p))
             {
               HD_p->addValue(1);
             }
          else
            {
              return;
            }
        }
      //HD_A -> nonHD_C:
      else if(variableA && C && !D && !variableD)
        {
          Voxel* moleculeC(reactvAC(variableA, C));
          if(moleculeC == NULL)
            {
              requeue();
              return;
            }
          else
            {
              variableA->addValue(-1);
              C->addMolecule(moleculeC);
            }
        }
      //HD_A -> nonHD_C + nonHD_D:
      else if(variableA && C && D)
        {
          Voxel* moleculeC(NULL);
          Voxel* moleculeD(NULL);
          Comp* compA(theSpatiocyteStepper->system2Comp(
                         variableA->getSuperSystem()));
          //Occupy C in a voxel of compartment C that adjoins compartment A
          //if A is a surface compartment:
          if(compA != C->getComp() && compA->dimension != 3)
            {
              moleculeC = C->getRandomAdjoiningCompVoxel(compA);
              if(moleculeC)
                {
                  moleculeD = D->getRandomAdjoiningVoxel(moleculeC, moleculeC);
                }
            }
          else if(compA != D->getComp() && compA->dimension != 3)
            {
              moleculeD = D->getRandomAdjoiningCompVoxel(compA);
              if(moleculeD)
                {
                  moleculeC = C->getRandomAdjoiningVoxel(moleculeD, moleculeD);
                }
            }
          else
            {
              moleculeC = C->getRandomCompVoxel();
              if(moleculeC)
                {
                  moleculeD = D->getRandomAdjoiningVoxel(moleculeC, moleculeC);
                }
            }
          if(moleculeC == NULL || moleculeD == NULL)
            {
              requeue();
              return;
            }
          variableA->addValue(-1);
          C->addMolecule(moleculeC);
          D->addMolecule(moleculeD);
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
          Voxel* molecule(reactvAC(variableA, nonHD_p));
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
  //theOrder = 2
  else
    {
      //HD + HD -> product(s)
      if(variableA && variableB)
        {
          //HD + HD -> HD: 
          if(variableC && !variableD && !D)
            {
              variableA->addValue(-1);
              variableB->addValue(-1);
              variableC->addValue(1);
            }
          //HD + HD -> nonHD: 
          else if(C && !variableD && !D)
            { 
              Voxel* moleculeC(reactvAvBC(C));
              if(moleculeC == NULL)
                {
                  requeue();
                  return;
                }
              variableA->addValue(-1);
              variableB->addValue(-1);
              C->addMolecule(moleculeC);
            }
          //HD + HD -> HD + HD: 
          else if(variableC && variableD)
            {
              variableA->addValue(-1);
              variableB->addValue(-1);
              variableC->addValue(1);
              variableD->addValue(1);
            }
          //HD + HD -> HD + nonHD: 
          //HD + HD -> nonHD + HD: 
          else if((variableC && D) || (C && variableD))
            {
              Variable* HD_p(variableC);
              Species* nonHD_p(D);
              if(variableD)
                {
                  HD_p = variableD;
                  nonHD_p = C;
                }
              Voxel* molecule(reactvAvBC(nonHD_p));
              if(molecule == NULL)
                {
                  requeue();
                  return;
                }
              variableA->addValue(-1);
              variableB->addValue(-1);
              nonHD_p->addMolecule(molecule);
              HD_p->addValue(1);
            }
          //HD + HD -> nonHD + nonHD: 
          else if(C && D)
            {
              Voxel* moleculeC(reactvAvBC(C));
              Voxel* moleculeD(NULL);
              if(moleculeC == NULL)
                {
                  moleculeD = reactvAvBC(D);
                  if(moleculeD)
                    {
                      moleculeC = C->getRandomAdjoiningVoxel(moleculeD,
                                                             moleculeD);
                    }
                }
              else
                { 
                  moleculeD = D->getRandomAdjoiningVoxel(moleculeC, moleculeC);
                }
              if(moleculeC == NULL || moleculeD == NULL)
                {
                  requeue();
                  return;
                }
              variableA->addValue(-1);
              variableB->addValue(-1);
              C->addMolecule(moleculeC);
              D->addMolecule(moleculeD);
            }
        }
      //HD + nonHD -> product(s)
      //nonHD + HD -> product(s)
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
              if(!reactACD(nonHD, C, D))
                {
                  return;
                }
              HD->addValue(-1);
            }
          //nonHD + HD -> nonHD:
          //HD + nonHD -> nonHD:
          else if(C && !D && !variableD)
            {
              if(!reactAC(nonHD, C))
                {
                  return;
                }
              HD->addValue(-1);
            }
          //HD + nonHD -> HD + nonHD:
          //HD + nonHD -> nonHD + HD:
          //nonHD + HD -> HD + nonHD:
          //nonHD + HD -> nonHD + HD:
          else if((variableC && D) || (C && variableD))
            {
              Variable* HD_p(variableC);
              Species* nonHD_p(D);
              if(variableD)
                {
                  HD_p = variableD;
                  nonHD_p = C;
                }
              if(!reactAC(nonHD, nonHD_p))
                {
                  return;
                }
              HD->addValue(-1);
              HD_p->addValue(1);
            }
        }
    }
  /*
  std::cout << "after" << std::endl;
  aValue1=theVariableReferenceVector[0].getVariable()->getValue();
  aValue2=theVariableReferenceVector[1].getVariable()->getValue();
  std::cout << "aValue1:" << aValue1 << " " << variableA->getValue();
  std::cout << " aValue2:" << aValue2 << " " << B->size() << std::endl;
  */
  ReactionProcess::fire();
}

//nonHD + nonHD -> nonHD
bool SpatiocyteNextReactionProcess::reactACD(Species* a, Species* c, Species* d)
{
  Voxel* moleculeA(a->getRandomMolecule());
  Voxel* moleculeC(NULL);
  Voxel* moleculeD(NULL);
  if(a->getVacantID() == c->getVacantID() || a->getID() == c->getVacantID())
    {
      moleculeC = moleculeA;
      moleculeD = d->getRandomAdjoiningVoxel(moleculeC, moleculeC);
      if(moleculeD == NULL)
        {
          requeue();
          return false;
        }
    }
  else if(a->getVacantID() == d->getVacantID() ||
          a->getID() == d->getVacantID())
    {
      moleculeD = moleculeA;
      moleculeC = c->getRandomAdjoiningVoxel(moleculeD, moleculeD);
      if(moleculeC == NULL)
        {
          requeue();
          return false;
        }
    }
  else
    {
      moleculeC = c->getRandomAdjoiningVoxel(moleculeA);
      if(moleculeC == NULL)
        {
          //Only proceed if we can find an adjoining vacant voxel
          //of nonND which can be occupied by C:
          requeue();
          return false;
        }
      moleculeD = d->getRandomAdjoiningVoxel(moleculeC, moleculeC);
      if(moleculeD == NULL)
        {
          requeue();
          return false;
        }
    }
  a->removeMolecule(moleculeA);
  c->addMolecule(moleculeC);
  d->addMolecule(moleculeD);
  return true;
}

//nonHD -> nonHD
bool SpatiocyteNextReactionProcess::reactAC(Species* a, Species* c)
{
  Voxel* moleculeA(a->getRandomMolecule());
  Voxel* moleculeC(NULL);
  if(a->getVacantID() == c->getVacantID() || a->getID() == c->getVacantID())
    {
      moleculeC = moleculeA;
    }
  else
    {
      moleculeC = c->getRandomAdjoiningVoxel(moleculeA);
      if(moleculeC == NULL)
        {
          //Only proceed if we can find an adjoining vacant voxel
          //of nonND which can be occupied by C:
          requeue();
          return false;
        }
    }
  a->removeMolecule(moleculeA);
  c->addMolecule(moleculeC);
  return true;
}

//nonHD (+Vacant[BindingSite]) -> nonHD
bool SpatiocyteNextReactionProcess::reactACbind(Species* a, Species* c)
{
  Voxel* moleculeA(a->getRandomMolecule());
  Voxel* moleculeC(NULL);
  moleculeC = c->getRandomAdjoiningVoxel(moleculeA, BindingSite);
  if(moleculeC == NULL)
    {
      //Only proceed if we can find an adjoining vacant voxel
      //of nonND which can be occupied by C:
      requeue();
      return false;
    }
  a->removeMolecule(moleculeA);
  c->addMolecule(moleculeC);
  return true;
}

//HD -> nonHD
Voxel* SpatiocyteNextReactionProcess::reactvAC(Variable* vA, Species* c)
{
  Voxel* moleculeC(NULL);
  Comp* compA(theSpatiocyteStepper->system2Comp(vA->getSuperSystem()));
  //Occupy C in a voxel of compartment C that adjoins compartment A
  //if A is a surface compartment:
  if(compA != c->getComp() && compA->dimension != 3)
    {
      moleculeC = c->getRandomAdjoiningCompVoxel(compA);
    }
  else
    {
      moleculeC = c->getRandomCompVoxel();
    }
  return moleculeC;
}

Comp* SpatiocyteNextReactionProcess::getComp2D(Species* c)
{
  Comp* compA(theSpatiocyteStepper->system2Comp(variableA->getSuperSystem()));
  Comp* compB(theSpatiocyteStepper->system2Comp(variableB->getSuperSystem()));
  Comp* comp2D(NULL);
  if(compA->dimension == 2)
    {
      comp2D = compA;
    }
  else if(compB->dimension == 2)
    {
      comp2D = compB;
    }
  //Occupy C in a voxel of compartment C that adjoins compartment A
  //if A is a surface compartment:
  if(comp2D != c->getComp() && comp2D != NULL)
    {
      return comp2D;
    }
  return NULL;
}

Voxel* SpatiocyteNextReactionProcess::reactvAvBC(Species* c)
{
  Voxel* moleculeC(NULL);
  Comp* aComp2D(getComp2D(c));
  if(aComp2D)
    {
      moleculeC = C->getRandomAdjoiningCompVoxel(aComp2D);
    }
  else
    {
      moleculeC = C->getRandomCompVoxel();
    }
  return moleculeC;
}

void SpatiocyteNextReactionProcess::initializeFourth()
{
  ReactionProcess::initializeFourth();
  if(p != -1)
    {
      return;
    }
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
  else if(variableB)
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
  double aVolume(0);
  double anArea(0);
  if(theOrder == 0)
    {
      double aSpace(0);
      if(SpaceC > 0)
        {
          aSpace = SpaceC;
          pFormula << "[aSpace:SpaceC:" << aSpace << "]";
        }
      else if(compC->dimension == 2)
        {
          aSpace = compC->actualArea;
          pFormula << "[aSpace:compC.Area:" << aSpace << "]";
        }
      else
        {
          aSpace = compC->actualVolume;
          pFormula << "[aSpace:compC.Volume:" << aSpace << "]";
        }
      p = k*aSpace;
      pFormula << "[k*aSpace:" << k << "*" << aSpace << "]";
    }
  else if(theOrder == 1) 
    {
      //Convert the unit m/s of k to 1/s for p if the reaction is a surface
      //adsorption reaction:
      if(compA->dimension == 3 && compC->dimension == 2)
        { 
          if(SpaceA > 0)
            {
              aVolume = SpaceA;
              pFormula << "[aVolume:SpaceA:" << aVolume << "]";
            }
          else
            {
              aVolume = compA->actualVolume;
              pFormula << "[aVolume:compA.Volume:" << aVolume << "]";
            }
          if(SpaceC > 0)
            {
              anArea = SpaceC;
              pFormula << "[anArea:SpaceC:" << anArea << "]";
            }
          else
            {
              anArea = compC->actualArea;
              pFormula << "[anArea:compC.Area:" << anArea << "]";
            }
          p = k*anArea/aVolume;
          pFormula << "[k*anArea/aVolume:" << k << "*" << anArea << "/"
            << aVolume << "]";
          return;
        }
      p = k;
      pFormula << "[k:" << k << "]";
    }
  else if(theOrder == 2)
    {
      //If there are two products that don't belong to the same compartment,
      //the reactants must also belong to different compartments:
      if((compD && compD != compC) && (compA == compB))
        {
          NEVER_GET_HERE;
        }
      //If volume + surface = k(volume)(surface) or
      //   volume + surface = k(surface)(volume) or
      //   surface + volume = k(volume)(surface) or
      //   surface + volume = k(surface)(volume)
      if((compD && (
        (compC->dimension == 3 && compD->dimension == 2 &&
         compA->dimension == 3 && compB->dimension == 2) ||
        (compC->dimension == 3 && compD->dimension == 2 &&
         compA->dimension == 2 && compB->dimension == 3) ||
        (compC->dimension == 2 && compD->dimension == 3 &&
         compA->dimension == 3 && compB->dimension == 2) ||
        (compC->dimension == 2 && compD->dimension == 3 &&
         compA->dimension == 2 && compB->dimension == 3))) ||
      //If volume (+volume) = k(volume)(volume) or
      //   surface (+surface) = k(volume)(surface) or
      //   surface (+surface) = k(surface)(volume)
         ((compC->dimension == 3 && compA->dimension == 3
          && compB->dimension == 3) ||
         (compC->dimension == 2 && compA->dimension == 3 
          && compB->dimension == 2) ||
         (compC->dimension == 2 && compA->dimension == 2 
          && compB->dimension == 3)))
        {
          if(compA->dimension == 3)
            {
              if(SpaceA > 0)
                {
                  aVolume = SpaceA;
                  pFormula << "[aVolume:SpaceA:" << aVolume << "]";
                }
              else
                {
                  aVolume = compA->actualVolume;
                  pFormula << "[aVolume:compA.Volume:" << aVolume << "]";
                }
            }
          else
            {
              if(SpaceB > 0)
                {
                  aVolume = SpaceB;
                  pFormula << "[aVolume:SpaceB:" << aVolume << "]";
                }
              else
                {
                  aVolume = compB->actualVolume;
                  pFormula << "[aVolume:compB.Volume:" << aVolume << "]";
                }
            }
          //unit of k is in m^3/s
          p = k/aVolume;
          pFormula << "[k/aVolume:" << k << "/" << aVolume << "]";
        }
      //If surface (+surface) = k(surface)(surface) or
      //   volume (+volume) = k(volume)(surface) or
      //   volume (+volume) = k(surface)(volume)
      else if((compC->dimension == 2 && compA->dimension == 2 
               && compB->dimension == 2) ||
              (compC->dimension == 3 && compA->dimension == 3 
               && compB->dimension == 2) ||
              (compC->dimension == 3 && compA->dimension == 2 
               && compB->dimension == 3))
        {
          if(compA->dimension == 2)
            {
              if(SpaceA > 0)
                {
                  anArea = SpaceA;
                  pFormula << "[anArea:SpaceA:" << anArea << "]";
                }
              else
                {
                  anArea = compA->actualArea;
                  pFormula << "[anArea:compA.Area:" << anArea << "]";
                }
            }
          else
            {
              if(SpaceB > 0)
                {
                  anArea = SpaceB;
                  pFormula << "[anArea:SpaceB:" << anArea << "]";
                }
              else
                {
                  anArea = compB->actualArea;
                  pFormula << "[anArea:compB.Area:" << anArea << "]";
                }
            }
          //unit of k is in m^2/s
          p = k/anArea;
          pFormula << "[k/anArea:" << k << "/" << anArea << "]";
        }
      else
        {
          NEVER_GET_HERE;
        }
      //A + A -> products
      if(getZeroVariableReferenceOffset() == 1)
        {
          p = k;
          pFormula << "[k:" << k << "]";
        }
    }
  else
    {
      NEVER_GET_HERE;
    } 
}

void SpatiocyteNextReactionProcess::printParameters()
{
  String aProcess(String(getPropertyInterface().getClassName()) + 
                                      "[" + getFullID().asString() + "]");
  std::cout << aProcess << std::endl;
  if(A)
    {
      std::cout << "  " << getIDString(A);
    }
  else if(variableA)
    {
      std::cout << "  " << getIDString(variableA);
    }
  if(B)
    {
      std::cout << " + " << getIDString(B);
    }
  else if(variableB)
    {
      std::cout << " + " << getIDString(variableB);
    }
  if(!A && !variableA)
    {
      if(C)
        {
          std::cout << "0 -> " << getIDString(C);
        }
      else if(variableC)
        {
          std::cout << "0 -> " << getIDString(variableC);
        }
    }
  else
    {
      if(C)
        {
          std::cout << " -> " << getIDString(C);
        }
      else if(variableC)
        {
          std::cout << " -> " << getIDString(variableC);
        }
    }
  if(D)
    {
      std::cout << " + " << getIDString(D);
    }
  else if(variableD)
    {
      std::cout << " + " << getIDString(variableD);
    }
  std::cout << " k:" << k << " p = " << pFormula.str() << " = " << p
    << " nextTime:" << getStepInterval() << std::endl;
}


GET_METHOD_DEF(Real, StepInterval, SpatiocyteNextReactionProcess)
{
  double step(getPropensity_R()*(-log(gsl_rng_uniform_pos(getStepper()->getRng()))));
  //std::cout << getFullID().asString() << " " << theTime <<  " next:" << theTime+step << " interval:" << step << std::endl; 
  return step;
}

