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
// based on GillespieProcess.hpp
// E-Cell Project, Institute for Advanced Biosciences, Keio University.
//


#ifndef __SpatiocyteNextReactionProcess_hpp
#define __SpatiocyteNextReactionProcess_hpp

#include <sstream>
#include <MethodProxy.hpp>
#include "ReactionProcess.hpp"

LIBECS_DM_CLASS(SpatiocyteNextReactionProcess, ReactionProcess)
{ 
  typedef MethodProxy<SpatiocyteNextReactionProcess, Real> RealMethodProxy; 
  typedef Real (SpatiocyteNextReactionProcess::*PDMethodPtr)(Variable*);
public:
  LIBECS_DM_OBJECT(SpatiocyteNextReactionProcess, Process)
    {
      INHERIT_PROPERTIES(ReactionProcess);
      PROPERTYSLOT_SET_GET(Real, SpaceA);
      PROPERTYSLOT_SET_GET(Real, SpaceB);
      PROPERTYSLOT_SET_GET(Real, SpaceC);
      PROPERTYSLOT_GET_NO_LOAD_SAVE(Real, Propensity);
    }
  SpatiocyteNextReactionProcess():
    initSizeA(0),
    initSizeB(0),
    initSizeC(0),
    initSizeD(0),
    SpaceA(0),
    SpaceB(0),
    SpaceC(0),
    theGetPropensityMethodPtr(RealMethodProxy::create<
            &SpatiocyteNextReactionProcess::getPropensity_ZerothOrder>()) {}
  virtual ~SpatiocyteNextReactionProcess() {}
  SIMPLE_SET_GET_METHOD(Real, SpaceA);
  SIMPLE_SET_GET_METHOD(Real, SpaceB);
  SIMPLE_SET_GET_METHOD(Real, SpaceC);
  virtual void initialize()
    {
      if(isInitialized)
        {
          return;
        }
      ReactionProcess::initialize();
      if(!(getOrder() == 0 || getOrder() == 1 || getOrder() == 2))
        {
          THROW_EXCEPTION(ValueError, 
                          String(getPropertyInterface().getClassName()) + 
                          "[" + getFullID().asString() + 
                          "]: Only zeroth, first or second order scheme " + 
                          "is allowed.");
        }
      if(variableA)
        {
          initSizeA = variableA->getValue();
        }
      if(variableB)
        {
          initSizeB = variableB->getValue();
        }
      if(variableC)
        {
          initSizeC = variableC->getValue();
        }
      if(variableD)
        {
          initSizeD = variableD->getValue();
        }
    }
  virtual void initializeSecond()
    {
      ReactionProcess::initializeSecond();
      if(variableA)
        {
          variableA->setValue(initSizeA);
        }
      if(variableB)
        {
          variableB->setValue(initSizeB);
        }
      if(variableC)
        {
          variableC->setValue(initSizeC);
        }
      if(variableD)
        {
          variableD->setValue(initSizeD);
        }
    }
  GET_METHOD(Real, Propensity)
    {
      Real aPropensity(theGetPropensityMethodPtr(this));
      if(aPropensity < 0.0)
        {
          THROW_EXCEPTION(SimulationError, "Variable value <= -1.0");
          return 0.0;
        }
      else
        {
          return aPropensity;
        }
    }
  GET_METHOD(Real, Propensity_R)
    {
      Real aPropensity(getPropensity());
      if(aPropensity > 0.0)
        {
          return 1.0/aPropensity;
        }
      else
        {
          return libecs::INF;
        }
    }
  virtual bool isContinuous() 
    {
      return true;
    }
  virtual GET_METHOD(Real, StepInterval);
  virtual void fire();
  virtual void initializeThird();
  virtual void printParameters();
protected:
  virtual void calculateOrder();
  Real getPropensity_ZerothOrder() 
    {
      return p;
    }
  Real getPropensity_FirstOrder() 
    {
      Real aValue(theVariableReferenceVector[0].getVariable()->getValue());
      if(aValue > 0.0)
        {
          return p*aValue;
        }
      else
        {
          return 0.0;
        }
    }
  Real getPropensity_SecondOrder_TwoSubstrates() 
    {
      Real aValue1(theVariableReferenceVector[0].getVariable()->getValue());
      Real aValue2(theVariableReferenceVector[1].getVariable()->getValue());
      if(aValue1 > 0.0 && aValue2 > 0.0)
        {
          return p*aValue1*aValue2;
        }
      else
        {
          return 0.0;
        }
    }
  Real getPropensity_SecondOrder_OneSubstrate() 
    {
      Real aValue(theVariableReferenceVector[0].getVariable()->getValue());
      //There must be two or more molecules:
      if(aValue > 1.0)
        {
          return p*aValue*(aValue-1.0);
        }
      else
        {
          return 0.0;
        }
    }
protected:
  double initSizeA;
  double initSizeB;
  double initSizeC;
  double initSizeD;
  double SpaceA;
  double SpaceB;
  double SpaceC;
  std::stringstream pFormula;
  RealMethodProxy theGetPropensityMethodPtr;  
};

inline void SpatiocyteNextReactionProcess::calculateOrder()
{
  ReactionProcess::calculateOrder();
  // set theGetPropensityMethodPtr
  if(getOrder() == 0) // no substrate
    {
      theGetPropensityMethodPtr = RealMethodProxy::create<
        &SpatiocyteNextReactionProcess::getPropensity_ZerothOrder>();
    }
  else if(getOrder() == 1)   // one substrate, first order.
    {
      theGetPropensityMethodPtr = RealMethodProxy::create<
        &SpatiocyteNextReactionProcess::getPropensity_FirstOrder>();
    }
  else if(getOrder() == 2)
    { 
      //Two unique substrate species, second order
      //A + B -> products:
      if(getZeroVariableReferenceOffset() == 2)
        {  
          theGetPropensityMethodPtr = RealMethodProxy::
            create<&SpatiocyteNextReactionProcess::
            getPropensity_SecondOrder_TwoSubstrates>();
        }
      //One substrate species, second order
      //A + A -> products:
      else
        {
          theGetPropensityMethodPtr = RealMethodProxy::
            create<&SpatiocyteNextReactionProcess::
            getPropensity_SecondOrder_OneSubstrate>();
        }
    }
  else
    {
      std::cout << "theOrder:" << getOrder() << std::endl;
      NEVER_GET_HERE;
    }
}

#endif /* __SpatiocyteNextReactionProcess_hpp */
