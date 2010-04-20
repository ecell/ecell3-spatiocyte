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
      PROPERTYSLOT_GET_NO_LOAD_SAVE(Real, Propensity);
    }
  SpatiocyteNextReactionProcess():
    theGetPropensityMethodPtr(RealMethodProxy::create<
                              &SpatiocyteNextReactionProcess::getZero>()),
    theGetMinValueMethodPtr(RealMethodProxy::create<
                            &SpatiocyteNextReactionProcess::getZero>()),
    theGetPDMethodPtr(&SpatiocyteNextReactionProcess::getPD_Zero) {}
  virtual ~SpatiocyteNextReactionProcess() {}
  virtual void initialize()
    {
      ReactionProcess::initialize();
      if(!(getOrder() == 1 || getOrder() == 2))
        {
          THROW_EXCEPTION(ValueError, 
                          String(getPropertyInterface().getClassName()) + 
                          "[" + getFullID().asString() + 
                          "]: Only first or second order scheme is allowed.");
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
  Real getPD(Variable* aVariable) 
    {
      return (this->*theGetPDMethodPtr)(aVariable);
    }
  virtual bool isContinuous() 
    {
      return true;
    }
  virtual GET_METHOD(Real, StepInterval);
  virtual void fire();
  virtual void initializeThird();
protected:
  virtual void calculateOrder();
  Real getZero() 
    {
      return 0.0;
    }
  Real getPD_Zero(Variable* aVariable) 
    {
      return 0.0;
    }
  Real getPropensity_FirstOrder() 
    {
      Real aValue(theVariableReferenceVector[0
                        ].getVariable()->getValue());
      if(aValue > 0.0)
        {
          return p*aValue;
        }
      else
        {
          return 0.0;
        }
    }
  Real getMinValue_FirstOrder() 
    {
      return theVariableReferenceVector[0].getVariable()->getValue();
    }
  Real getPD_FirstOrder(Variable* aVariable) 
    {
      if(theVariableReferenceVector[0].getVariable() == aVariable)
        {
          return p;
        }
      else
        {
          return 0.0;
        }
    }
  Real getPropensity_SecondOrder_TwoSubstrates() 
    {
      Real aValue1(theVariableReferenceVector[0
                         ].getVariable()->getValue());
      Real aValue2(theVariableReferenceVector[1
                         ].getVariable()->getValue());
      if(aValue1 > 0.0 && aValue2 > 0.0)
        {
          return p*aValue1*aValue2;
        }
      else
        {
          return 0.0;
        }
    }
  Real getMinValue_SecondOrder_TwoSubstrates() 
    {
      Real aFirstValue(theVariableReferenceVector[0
                             ].getVariable()->getValue());
      Real aSecondValue(theVariableReferenceVector[1
                              ].getVariable()->getValue());
      return fmin(aFirstValue, aSecondValue);
    }
  Real getPD_SecondOrder_TwoSubstrates(Variable* aVariable) 
    {
      if(theVariableReferenceVector[0].getVariable() == aVariable)
        {
          return p*theVariableReferenceVector[1].getVariable()->getValue();
        }
      else if(theVariableReferenceVector[1].getVariable() == aVariable)
        {
          return p*theVariableReferenceVector[0].getVariable()->getValue();
        }
      else
        {
          return 0.0;
        }
    }
  Real getPropensity_SecondOrder_OneSubstrate() 
    {
      Real aValue(theVariableReferenceVector[0
                        ].getVariable()->getValue());
      if (aValue > 1.0) // there must be two or more molecules
        {
          return p*0.5*aValue*(aValue-1.0);
        }
      else
        {
          return 0.0;
        }
    }
  Real getMinValue_SecondOrder_OneSubstrate() 
    {
      return theVariableReferenceVector[0].getVariable()->getValue()*0.5;
    }
  Real getPD_SecondOrder_OneSubstrate(Variable* aVariable) 
    {
      if(theVariableReferenceVector[0].getVariable() == aVariable)
        {
          Real aValue(theVariableReferenceVector[0
                            ].getVariable()->getValue());
          if(aValue > 0.0) // there must be at least one molecule
            {
              return p*0.5*(2.0*aValue-1.0);
            }
          else
            {
              return 0.0;
            }
        }
      else
        {
          return 0.0;
        }
    }
protected:
  RealMethodProxy theGetPropensityMethodPtr;  
  RealMethodProxy theGetMinValueMethodPtr;
  PDMethodPtr     theGetPDMethodPtr;
};

inline void SpatiocyteNextReactionProcess::calculateOrder()
{
  ReactionProcess::calculateOrder();
  for(VariableReferenceVector::iterator
      i(theVariableReferenceVector.begin());
      i != theVariableReferenceVector.end(); ++i)
    {
      VariableReference aVariableReference(*i);
      Integer aCoefficient( aVariableReference.getCoefficient() );
      // here assume aCoefficient != 0
      if(aCoefficient == 0)
        {
          THROW_EXCEPTION(InitializationFailed,
                          "[" + getFullID().asString() + 
                          "]: Zero stoichiometry is not allowed.");
        }
    }
  // set theGetPropensityMethodPtr and theGetMinValueMethodPtr
  if(getOrder() == 0)   // no substrate
    {
      theGetPropensityMethodPtr =
        RealMethodProxy::create<&SpatiocyteNextReactionProcess::getZero>();
      theGetMinValueMethodPtr =
        RealMethodProxy::create<&SpatiocyteNextReactionProcess::getZero>();
      theGetPDMethodPtr = &SpatiocyteNextReactionProcess::getPD_Zero;
    }
  else if(getOrder() == 1)   // one substrate, first order.
    {
      theGetPropensityMethodPtr =
        RealMethodProxy::create<
        &SpatiocyteNextReactionProcess::getPropensity_FirstOrder>();
      theGetMinValueMethodPtr =
        RealMethodProxy::create<
        &SpatiocyteNextReactionProcess::getMinValue_FirstOrder>();
      theGetPDMethodPtr = &SpatiocyteNextReactionProcess::getPD_FirstOrder;
    }
  else if(getOrder() == 2)
    {
      if(getZeroVariableReferenceOffset() == 2) // 2 substrates, 2nd order
        {  
          theGetPropensityMethodPtr = RealMethodProxy::
            create<&SpatiocyteNextReactionProcess::
            getPropensity_SecondOrder_TwoSubstrates>();
          theGetMinValueMethodPtr = RealMethodProxy::create<
            &SpatiocyteNextReactionProcess::
            getMinValue_SecondOrder_TwoSubstrates>();
          theGetPDMethodPtr = 
            &SpatiocyteNextReactionProcess::getPD_SecondOrder_TwoSubstrates;
        }
      else // one substrate, second order (coeff == -2)
        {
          theGetPropensityMethodPtr = RealMethodProxy::
            create<&SpatiocyteNextReactionProcess::
            getPropensity_SecondOrder_OneSubstrate>();
          theGetMinValueMethodPtr = RealMethodProxy::create<
            &SpatiocyteNextReactionProcess::
            getMinValue_SecondOrder_OneSubstrate>();
          theGetPDMethodPtr =
            &SpatiocyteNextReactionProcess::getPD_SecondOrder_OneSubstrate;
        }
    }
  else
    {
      //FIXME: generic functions should come here.
      theGetPropensityMethodPtr =
        RealMethodProxy::create<&SpatiocyteNextReactionProcess::getZero>();
      theGetPropensityMethodPtr =
        RealMethodProxy::create<&SpatiocyteNextReactionProcess::getZero>();
      theGetPDMethodPtr = &SpatiocyteNextReactionProcess::getPD_Zero;
    }
}

#endif /* __SpatiocyteNextReactionProcess_hpp */
