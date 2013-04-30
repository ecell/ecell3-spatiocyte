import numpy as np
from scipy import *
from scipy.optimize import curve_fit


def fit_linear(parameters, values):
	a, b = parameters
	x = values
	return (a + b * x)

def fit_cuadratic(parameters, values):
	x = values
	x0, a, b, c = parameters
	return(a *(x-x0)**2 + b*(x-x0) + c )

def fit_gauss(parameters, values):
	a, x0, w, y0 = parameters
	x = values
	return(a * np.exp(-2*(x-x0)**2/w**2) + y0)

def fit_lorentz(parameters, values):
	a, x0, w = parameters
	x = values
	return( a * w **2 / ((x-x0)**2 + w**2) )

def fit_boltzman(parameters, values):
	x0, a1, a2, dx = parameters
	x = values
	return((a1 - a2)/(1 + np.exp((x - x0)/dx)) + a2)

def fit_logistic(parameters, values):
	x0, a1, a2, p = parameters
	x = values
	return((a1 - a2)/(1 + (x/x0)**p) + a2)

def fit_expgrow(parameters, values):
	x0, y0, a, t = parameters
	x = values
	return(y0 + a * np.exp((x - x0)/t))

def fit_expassoc(parameters, values):
	y0, a1, t1, a2, t2 = parameters
	x = values
	return(y0 + a1 * (1 + np.exp(-x/t1)) + a2 * (1 + np.exp(-x/t2)))

def fit_hyperbl(parameters, values):
	p1, p2 = parameters
	x = values
	return(p1 * x/ (p2 + x))

def fit_pulse(parameters, values):
	x0, y0, a, t1, t2 = parameters
	x = values
	return(y0 + a * (1 + np.exp(-(x - x0)/t1)) * np.exp(-(x - x0)/t2))

def fit_rational0(parameters, values):
	a, b, c = parameters
	x = values
	return((b + c*x)/(1 + a*x))

def fit_sine(parameters, values):
	x0, a, w = parameters
	x = values
	return(a * sin(pi*(x - x0)/w))

def fit_gaussamp(parameters, values):
	x = values
	a, x0, w, y0 = parameters
	return( a * np.exp(-(x - x0)**2/(2*w**2)) + y0)
#	a0, a1, a2, a3, a4, a5 = parameters
#	return(a0 * exp(-(x-a1)**2/a2**2/2) *( a3*x**2 + a4 * x + a5 ))

def fit_allometric(parameters, values):
	a, x0, c = parameters
	x = values
	return(a * (x-x0)**2 + c)

def fit_expdecay(x, a, t):
  x0 = 0
  y0 = 3.01
  return(y0 + a*np.exp(-(x-x0)/t))

def func(x, a, b, c):
  return a*np.exp(-b*x) + c

x = array([0,10,20,30, 50])
print "x:", x
yn = array([3.01, 1.42, 1.16, 0.74, 0.63])
print "yin:", yn
popt, pcov = curve_fit(func, x, yn)

print "a b c:", popt

y = func(x, popt[0], popt[1], popt[2])
print "y:", y
