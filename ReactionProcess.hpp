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


#ifndef __ReactionProcess_hpp
#define __ReactionProcess_hpp

#include "SpatiocyteProcess.hpp"
#include "ReactionProcessInterface.hpp"

LIBECS_DM_CLASS_EXTRA_1(ReactionProcess, SpatiocyteProcess, ReactionProcessInterface)
{ 
public:
  LIBECS_DM_OBJECT(ReactionProcess, Process)
    {
      INHERIT_PROPERTIES(SpatiocyteProcess);
      PROPERTYSLOT_SET_GET(Real, k);
      PROPERTYSLOT_SET_GET(Real, p);
      PROPERTYSLOT_SET_GET(Integer, SearchVacant);
      PROPERTYSLOT_GET_NO_LOAD_SAVE(Integer, Order);
    }
  ReactionProcess():
    SearchVacant(-1),
    theOrder(0),
    k(-1),
    p(-1),
    A(NULL),
    B(NULL),
    C(NULL),
    D(NULL),
    E(NULL),
    variableA(NULL),
    variableB(NULL),
    variableC(NULL),
    variableD(NULL), 
    variableE(NULL),
    moleculeA(NULL),
    moleculeB(NULL),
    moleculeC(NULL),
    moleculeD(NULL),
    moleculeE(NULL),
    moleculeP(NULL),
    moleculeS(NULL) {}
  virtual ~ReactionProcess() {}
  SIMPLE_SET_GET_METHOD(Real, k);
  SIMPLE_SET_GET_METHOD(Real, p);
  SIMPLE_SET_GET_METHOD(Integer, SearchVacant);
  virtual bool isInterrupting(Process*);
  virtual void fire()
    {
      const Time aCurrentTime(theTime); // do this only for the Processes in Q
      requeue(); //theTop in thePriorityQueue is still this process since
      //we have not interrupted other processes to update their queue. 
      //So it is valid to call requeue, which only requeues theTop process, 
      //assuming it to be this process.
      for(std::vector<SpatiocyteProcess*>::const_iterator 
          i(theInterruptingProcesses.begin());
          i!=theInterruptingProcesses.end(); ++i)
        {
          (*i)->substrateValueChanged(aCurrentTime);
        }
    }
  GET_METHOD(Integer, Order)
    {
      return theOrder;
    }
  virtual void initialize()
    {
      if(isInitialized)
        {
          return;
        }
      SpatiocyteProcess::initialize();
      //SearchVacant of each reaction is set according to the master
      //SearchVacant property of the SpatiocyteStepper if it is not
      //specified at the reaction level:
      if(SearchVacant == -1)
        {
          SearchVacant = theSpatiocyteStepper->getSearchVacant();
        }
      declareUnidirectional();
      calculateOrder();
    }
  virtual void initializeSecond()
    {
      SpatiocyteProcess::initializeSecond();
      theInterruptingProcesses.resize(0);
    }
  virtual void initializeFourth()
    {
      SpatiocyteProcess::initializeFourth();
    }
  //This method is called whenever the substrate (-1 and 0 coefficients)
  //value of the process listed in isInterrupting method has changed:
  virtual void setInterrupt(std::vector<Process*> const &aProcessList,
                            Process* aProcess)
    {
      for(std::vector<Process*>::const_iterator i(aProcessList.begin());
          i != aProcessList.end(); ++i)
        {
          if(aProcess != (*i) && isInterrupting(*i))
            {
              theInterruptingProcesses.push_back(
                                     dynamic_cast<SpatiocyteProcess*>(*i));

            }
        }
    }
  virtual Species* getA()
    {
      return A;
    }
  virtual Species* getB()
    {
      return B;
    }
  virtual Species* getC()
    {
      return C;
    }
  virtual Species* getD()
    {
      return D;
    }
  virtual Species* getE()
    {
      return E;
    }
  virtual Voxel* getMoleculeA()
    {
      return moleculeA;
    }
  virtual Voxel* getMoleculeB()
    {
      return moleculeB;
    }
  virtual Voxel* getMoleculeC()
    {
      return moleculeC;
    }
  virtual Voxel* getMoleculeD()
    {
      return moleculeD;
    }
  virtual Voxel* getMoleculeE()
    {
      return moleculeE;
    }
  virtual Voxel* getMoleculeP()
    {
      return moleculeP;
    }
  virtual Voxel* getMoleculeS()
    {
      return moleculeS;
    }
  virtual void addSubstrateInterrupt(Species* aSpecies, Voxel* aMolecule) {}
  virtual void removeSubstrateInterrupt(Species* aSpecies, Voxel* aMolecule) {}
protected:
  virtual void calculateOrder();
protected:
  int SearchVacant;
  int theOrder;
  double k;
  double p;
  //Species are for non HD species:
  Species* A;
  Species* B;
  Species* C;
  Species* D;
  Species* E;
  //Variables are for HD species:
  Variable* variableA;
  Variable* variableB;
  Variable* variableC;
  Variable* variableD;
  Variable* variableE;
  Voxel* moleculeA;
  Voxel* moleculeB;
  Voxel* moleculeC;
  Voxel* moleculeD;
  Voxel* moleculeE;
  Voxel* moleculeP;
  Voxel* moleculeS;
  std::vector<SpatiocyteProcess*> theInterruptingProcesses;
};

inline void ReactionProcess::calculateOrder()
{ 
  theOrder = 0;
  for(VariableReferenceVector::iterator 
      i(theSortedVariableReferences.begin());
      i != theSortedVariableReferences.end(); ++i)
    {
      const int aCoefficient((*i).getCoefficient());
      Variable* aVariable((*i).getVariable());
      if(aCoefficient < 0)
        {
          theOrder -= aCoefficient; 
          //The first reactant, A:
          if(A == NULL && variableA == NULL)
            {
              if(aVariable->getName() == "HD")
                {
                  variableA = aVariable;
                }
              else
                {
                  A = theSpatiocyteStepper->getSpecies(aVariable);
                }
            }
          //The second reactant, B:
          else
            {
              if(aVariable->getName() == "HD")
                {
                  variableB = aVariable;
                }
              else
                {
                  B = theSpatiocyteStepper->getSpecies(aVariable);
                }
            }
        }
      else if(aCoefficient > 0)
        {
          //The first product, C:
          if(C == NULL && variableC == NULL)
            {
              if(aVariable->getName() == "HD")
                {
                  variableC = aVariable;
                }
              else
                {
                  C = theSpatiocyteStepper->getSpecies(aVariable);
                }
            }
          //The second product, D:
          else
            {
              if(aVariable->getName() == "HD")
                {
                  variableD = aVariable;
                }
              else
                {
                  D = theSpatiocyteStepper->getSpecies(aVariable);
                }
            }
        }
      //aCoefficient == 0:
      else
        {
          if(aVariable->getName() == "HD")
            {
              variableE = aVariable;
            }
          else
            {
              E = theSpatiocyteStepper->getSpecies(aVariable);
            }
        }
    }
} 

bool ReactionProcess::isInterrupting(Process* aProcess)
{
  //List all the processes here that need to be notified when their
  //substrateValueChanged:
  if(aProcess->getPropertyInterface().getClassName() ==
     "SpatiocyteNextReactionProcess" /* ||
     aProcess->getPropertyInterface().getClassName() == "DiffusionProcess"*/) 
    {
      //First get the unique variable pointers of this process:
      std::vector<Variable*> aVariableList;
      for(VariableReferenceVector::iterator
          i(theVariableReferenceVector.begin());
          i != theVariableReferenceVector.end(); ++i)
        {
          std::vector<Variable*>::const_iterator j(aVariableList.begin());
          while(j!=aVariableList.end())
            {
              if((*i).getVariable() == (*j))
                {
                  break;
                }
              ++j;
            }
          if(j == aVariableList.end())
            {
              aVariableList.push_back((*i).getVariable());
            }
        }
      //Find out if the values of the unique variables will be changed
      //by this process, i.e, netCoefficient != 0:
      std::vector<int> aNetCoefficientList;
      aNetCoefficientList.resize(aVariableList.size());
      for(std::vector<int>::iterator i(aNetCoefficientList.begin());
          i!=aNetCoefficientList.end(); ++i)
        {
          (*i) = 0;
        }
      for(VariableReferenceVector::iterator
          i(theVariableReferenceVector.begin());
          i != theVariableReferenceVector.end(); ++i)
        {
          for(std::vector<Variable*>::const_iterator j(aVariableList.begin());
              j!=aVariableList.end(); ++j)
            {
              if((*i).getVariable() == (*j))
                {
                  aNetCoefficientList[j-aVariableList.begin()] +=
                    (*i).getCoefficient();
                }
            }
        }
      //Check if any variable with netCoefficient != 0 is a substrate
      //of aProcess:
      VariableReferenceVector
        aVariableReferenceVector(aProcess->getVariableReferenceVector()); 
      for(VariableReferenceVector::iterator
          i(aVariableReferenceVector.begin());
          i != aVariableReferenceVector.end(); ++i)
        {
          if((*i).isAccessor())
            {
              for(std::vector<Variable*>::const_iterator j(aVariableList.begin());
                  j!=aVariableList.end(); ++j)
                {
                  if((*i).getVariable() == (*j) && 
                     aNetCoefficientList[j-aVariableList.begin()])
                    {
                      return true;
                    }
                }
            }
        }
    }
  return false;
}


#endif /* __ReactionProcess_hpp */
