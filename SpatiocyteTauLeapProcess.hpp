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


#ifndef __SpatiocyteTauLeapProcess_hpp
#define __SpatiocyteTauLeapProcess_hpp

#include <SpatiocyteNextReactionProcess.hpp> 


LIBECS_DM_CLASS(SpatiocyteTauLeapProcess, SpatiocyteNextReactionProcess)
{ 
  typedef double (SpatiocyteTauLeapProcess::*RealMethod)(double);

public:
  LIBECS_DM_OBJECT(SpatiocyteTauLeapProcess, Process)
    {
      INHERIT_PROPERTIES(SpatiocyteNextReactionProcess);
    }
  SpatiocyteTauLeapProcess():
    isParent(false),
    n(10),
    n_c(10),
    epsilon(0.03) {}
  virtual ~SpatiocyteTauLeapProcess() {}
  virtual void initialize()
    {
      SpatiocyteNextReactionProcess::initialize();
      const std::vector<Process*>& aProcesses(
                      theSpatiocyteStepper->getProcessVector());
      for(unsigned i(0); i != aProcesses.size(); ++i)
        {
          SpatiocyteTauLeapProcess* aProcess(
               dynamic_cast<SpatiocyteTauLeapProcess*>(aProcesses[i]));
          if(aProcess)
            {
              if(aProcess->getIsParent())
                {
                  R.resize(0);
                  return;
                }
              R.push_back(aProcess);
            }
        }
      isParent = true;
      initPoisson();
    }
  virtual void initializeFirst()
    {
      SpatiocyteNextReactionProcess::initializeFirst();
      if(isParent)
        {
          for(unsigned j(0); j != R.size(); ++j)
            {
              if(R[j] != this)
                {
                  if(R[j]->getIsExternInterrupted())
                    {
                      isExternInterrupted = true;
                      R[j]->setExternInterrupted(false);
                    }
                  if(R[j]->getN_c() != 10)
                    {
                      n_c = R[j]->getN_c();
                    }
                  if(R[j]->getEpsilon() != 0.03)
                    {
                      epsilon = R[j]->getEpsilon();
                    }
                }
              R[j]->setS_Index(S_rs, HOR, g);
            }
        }
      else
        {
          isPriorityQueued = false;
        }
    }
  bool getIsParent()
    {
      return isParent;
    }
  unsigned getL()
    {
      unsigned L(UINT_MAX);
      for(unsigned i(0); i != S_netNeg.size(); ++i)
        {
          const unsigned val(S_netNeg[i]->getValue()/(-v_netNeg[i]));
          if(val < L)
            {
              L = val;
            }
        }
      return L;
    }
  virtual double getInterval(double aCurrentTime)
    {
      return getNewInterval(aCurrentTime);
    }
  virtual double getNewInterval(double aCurrentTime)
    {
      double a0(0);
      double a0_c(0);
      tau1 = getTau(a0, a0_c);
      if(a0)
        {
          if(tau1 < n/a0)
            {
              std::cout << "do SSA" << std::endl;
              theState = 0;
              return libecs::INF;
            }
          tau2 = (a0_c == 0)? libecs::INF : -log(theRng->FixedU())/a0_c;
          if(tau1 < tau2)
            {
              theState = 1;
              return tau1;
            }
          theState = 2;
          return tau2;
        }
      return libecs::INF;
    }
  virtual void fire()
    {
      switch(theState)
        {
        case 0:
          std::cout << "firing SSA" << std::endl;
          break;
        case 1:
          fireNonCritical(tau1);
          break;
        case 2:
          fireCritical();
          fireNonCritical(tau2);
          break;
        }
      SpatiocyteNextReactionProcess::requeue();
    }
  void fireChild()
    {
      SpatiocyteNextReactionProcess::fire();
    }
  void fireNonCritical(const double aTau)
    {
      for(unsigned j(0); j != R.size(); ++j)
        {
          if(!R[j]->getIsCritical())
            {
              const unsigned K(poisson(R[j]->getPropensity()*aTau));
              for(unsigned i(0); i != K; ++i)
                {
                  R[j]->fireChild();
                }
            }
        }
    }
  void fireCritical()
    {
      unsigned j(theRng->Integer(R.size()));
      for(; j != R.size(); ++j)
        {
          if(R[j]->getIsCritical())
            {
              R[j]->fireChild();
              return;
            }
        }
      for(unsigned i(0); i != j; ++i)
        {
          if(R[i]->getIsCritical())
            {
              R[i]->fireChild();
              return;
            }
        }
    }
  void initPoisson()
    {
      poisson_table = {
        0.0,
        0.0,
        0.69314718055994529,
        1.7917594692280550,
        3.1780538303479458,
        4.7874917427820458,
        6.5792512120101012,
        8.5251613610654147,
        10.604602902745251,
        12.801827480081469};
    }
  unsigned poisson(double _mean)
    {
      using std::floor;
      using std::abs;
      using std::log; 
      using std::sqrt;
      _ptrd.smu = sqrt(_mean);
      _ptrd.b = 0.931 + 2.53 * _ptrd.smu;
      _ptrd.a = -0.059 + 0.02483 * _ptrd.b;
      _ptrd.inv_alpha = 1.1239 + 1.1328 / (_ptrd.b - 3.4);
      _ptrd.v_r = 0.9277 - 3.6224 / (_ptrd.b - 2);
      while(true)
        {
          double u;
          double v(theRng->Fixed());
          if(v <= 0.86*_ptrd.v_r)
            {
              u = v/_ptrd.v_r - 0.43;
              return static_cast<unsigned>(floor(
                (2*_ptrd.a/(0.5-abs(u)) + _ptrd.b)*u + _mean + 0.445));
            } 
          if(v >= _ptrd.v_r)
            {
              u = theRng->Fixed() - 0.5;
            }
          else
            {
              u = v/_ptrd.v_r - 0.93;
              u = ((u < 0)? -0.5 : 0.5) - u;
              v = theRng->Fixed()*_ptrd.v_r;
            } 
          const double us(0.5 - abs(u));
          if(us < 0.013 && v > us)
            {
              continue;
            } 
          const double K(floor((2*_ptrd.a/us + _ptrd.b)*u+_mean+0.445));
          v = v*_ptrd.inv_alpha/(_ptrd.a/(us*us) + _ptrd.b); 
          const double log_sqrt_2pi(0.91893853320467267); 
          if(K >= 10)
            {
              if(log(v*_ptrd.smu) <= (K + 0.5)*log(_mean/K)
                               - _mean
                               - log_sqrt_2pi
                               + K
                               - (1/12. - (1/360. - 1/(1260.*K*K))/(K*K))/K)
                {
                  return static_cast<unsigned>(K);
                }
            }
          else if(K >= 0)
            {
              if(log(v) <= K*log(_mean)
                           - _mean
                           - poisson_table[static_cast<unsigned>(K)])
                {
                  return static_cast<unsigned>(K);
                }
            }
        }
    }
  virtual void requeue() {}
  double getTau(double& a0, double& a0_c)
    {
      double tau(libecs::INF);
      std::vector<double> mu(S_rs.size(), 0);
      std::vector<double> sigma(S_rs.size(), 0);
      for(unsigned j(0); j != R.size(); ++j)
        {
          const double a(R[j]->getPropensity());
          if(a)
            {
              a0 += a;
              //If the reaction is non-critical:
              if(R[j]->getL() >= n_c)
                {
                  R[j]->addMuSigma(mu, sigma, a);
                }
              else
                {
                  a0_c += a;
                  R[j]->setIsCritical();
                }
            }
        }
      for(unsigned i(0); i != S_rs.size(); ++i)
        {
          const double x(S_rs[i]->getValue());
          const double ex_g(std::max(epsilon*x/(this->*g[i])(x), 1.0));
          const double tmp(std::min(ex_g/fabs(mu[i]), ex_g*ex_g/sigma[i]));
          if(tmp < tau)
            {
              tau = tmp;
            }
        }
      return tau;
    }
  void setIsCritical()
    {
      isCritical = true;
    }
  bool getIsCritical()
    {
      return isCritical;
    }
  void addMuSigma(std::vector<double>& mu_p, std::vector<double>& sigma_p,
                  const double& a)
    {
      isCritical = false;
      for(unsigned i(0); i != S_netNeg.size(); ++i)
        {
          const double aMu(v_netNeg[i]*a);
          mu_p[S_index[i]] += aMu;
          sigma_p[S_index[i]] += v_netNeg[i]*aMu;
        }
    }
  double getPropensity()
    {
      return (this->*thePropensityMethod)();
    }
  double getEpsilon()
    {
      return epsilon;
    }
  unsigned getN_c()
    {
      return n_c;
    }
  void setExternInterrupted(bool value)
    {
      isExternInterrupted = value;
    }
  void setS_Index(std::vector<Variable*>& aS, 
                  std::vector<unsigned>& aHOR,
                  std::vector<RealMethod>& aG)
    {
      setNetCoefficients();
      S_index.resize(S_netNeg.size());
      for(unsigned i(0); i != S_netNeg.size(); ++i)
        {
          Variable* aVariable(S_netNeg[i]);
          std::vector<Variable*>::iterator iter(std::find(aS.begin(), aS.end(),
                                                          aVariable));
          if(iter != aS.end())
            {
              const unsigned index(iter-aS.begin());
              S_index[i] = index;
              set_g(aHOR[index], aG[index], v_neg[i]);
            }
          else
            {
              S_index[i] = aS.size();
              aS.push_back(aVariable);
              aHOR.push_back(0);
              aG.push_back(NULL);
              set_g(aHOR.back(), aG.back(), v_neg[i]);
            }
        }
    }
  double g_order_1(double x)
    {
      return 1;
    }
  double g_order_2_1(double x)
    {
      return 2;
    }
  double g_order_2_2(double x)
    {
      return 2+1/(x-1);
    }
  double g_order_3_1(double x)
    {
      return 3;
    }
  double g_order_3_2(double x)
    {
      return 3/2*(2+1/(x-1));
    }
  double g_order_3_3(double x)
    {
      return 3+1/(x-1)+2/(x-2);
    }
  void setNetCoefficients()
    {
      std::vector<int> v;
      //First get the unique Variables of this process, S_net:
      for(VariableReferenceVector::const_iterator 
          i(theVariableReferenceVector.begin());
          i != theVariableReferenceVector.end(); ++i)
        {
          Variable* aVariable((*i).getVariable());
          std::vector<Variable*>::iterator iter(std::find(S_net.begin(),
                          S_net.end(), aVariable));
         if(iter == S_net.end())
            {
              S_net.push_back(aVariable);
              if((*i).getCoefficient() < 0)
                {
                  v.push_back((*i).getCoefficient());
                }
              else
                {
                  v.push_back(0);
                }
            }
          else if((*i).getCoefficient() < 0)
            {
              v[iter-S_net.begin()] += (*i).getCoefficient();
            }
        }
      //Find out if the values of the unique variables will be changed
      //by the Process aProcess, i.e., netCoefficient != 0:
      v_net.resize(S_net.size(), 0);
      for(VariableReferenceVector::const_iterator
          i(theVariableReferenceVector.begin());
          i != theVariableReferenceVector.end(); ++i)
        {
          for(std::vector<Variable*>::const_iterator
              j(S_net.begin()); j != S_net.end(); ++j)
            {
              if((*i).getVariable() == (*j))
                {
                  v_net[j-S_net.begin()] += (*i).getCoefficient();
                }
            }
        }
      for(unsigned i(0); i != v_net.size(); ++i)
        {
          if(v_net[i] < 0)
            {
              v_netNeg.push_back(v_net[i]);
              S_netNeg.push_back(S_net[i]);
              v_neg.push_back(v[i]);
            }
        }
    }
  void set_g(unsigned&, RealMethod&, const int);
private:
  bool isCritical;
  bool isParent;
  unsigned n;
  unsigned n_c;
  unsigned theState;
  double epsilon;
  double tau1;
  double tau2;
  std::vector<unsigned> S_index;
  std::vector<unsigned> HOR;
  std::vector<int> v_neg;
  std::vector<int> v_net;
  std::vector<int> v_netNeg;
  std::vector<RealMethod> g; 
  std::vector<SpatiocyteTauLeapProcess*> R;
  std::vector<Variable*> S_net;
  std::vector<Variable*> S_netNeg;
  std::vector<Variable*> S_rs;
  std::vector<double> poisson_table;
  struct
    {
      double v_r;
      double a;
      double b;
      double smu;
      double inv_alpha;
    } _ptrd;
};

#endif /* __SpatiocyteTauLeapProcess_hpp */

