'''
Check diagrammatic Monte Carlo test for
self-energy up to the second order

Last updates: April-12-2019
'''

###################################
##########Import libraries#########
###################################
#Begin
from math import *
from random import random
from random import randint
import numpy as np
import scipy.integrate 

#End
###################################

###################################
###########Fermi Function##########
###################################
#Begin
def fermi(e):  # e: dispersion energy, B: inverse temperature
  val = 1.0/(exp(beta*e)+1.0) 
  return val
#End
###################################

###################################
###########Bose Function##########
###################################
#Begin
def bose(e):  # e: dispersion energy, B: inverse temperature
  if e!=0:
    val = 1.0/(exp(beta*(e))-1.0) 
  else:          # if e==0 then use a regulator
    reg=0.01
    val = 1.0/(exp(beta*e-reg)-1.0) 
  return val
#End
###################################

###################################
############Dispersion#############
###################################
#Begin
def E(p): # p: momentum, m: chemical potential, t: hopping amplitude, a: lattice constant
  val = -2*t*(cos(p[0]*lat_const)+cos(p[1]*lat_const))-mu
  return val
#End
###################################

###################################
##########Integrand zeroth#########
###################################
#Begin
def Integrand_zeroth_order(p):
  return 1
#End
###################################

###################################
########Integrand Hartree##########
###################################
#Begin
def Integrand_1st_order(p): # p: momentum, m: chemical potential, t: hopping amplitude, a: lattice constant, B: inverse temperature, u: Hubbard potential
  dis_en=E(p)
  val=U*fermi(dis_en)/2/2/pi/pi
  return val
###################################

############Integrand##############
###################################
#Begin
def Integrand_2nd_order(p): # p: momentum, m: chemical potential, t: hopping amplitude, a: lattice constant, B: inverse temperature, u: Hubbard potential
  k=[p[0],p[1]]
  e1=E(k)
  k = [p[2],p[3]]
  e2=E(k)
  k=[px_ext+p[2]-p[0],py_ext+p[3]-p[1]]
  e3=E(k)
  val=fermi(e1)-fermi(-e3)
  val=val*(fermi(e2)+bose(e1+e3)) 
  val=-U*U*val/2/2/2/2/pi/pi/pi/pi  # Only the prefactor
  val=val/(1j*nu_ext-e1-e3+e2)
  return val
#End
###################################

###################################
##########Integrand real###########
###################################
#Begin
def integrand_real(p):
  ORDER=len(p)/2
  if ORDER==0:
    return Integrand_zeroth_order(p)
  if ORDER==1:
    return Integrand_1st_order(p)
  if ORDER==2:
    return Integrand_2nd_order(p).real
#End
###################################

###################################
##########Integrand imag###########
###################################
#Begin
def integrand_imag(p):
  ORDER=len(p)/2
  if ORDER==0:
    return Integrand_zeroth_order(p)
  if ORDER==1:
    return Integrand_1st_order(p)
  if ORDER==2:
    return Integrand_2nd_order(p).imag
#End
###################################

###################################
######Read external variables######
###################################
#Begin
data = np.loadtxt('ext_vars.dat')  # reading external variables from a file
count = len(open('ext_vars.dat').readlines())
if count==1:
  EXT_VARS = [None]*count
  EXT_VARS[0]=data
if count>1:
  EXT_VARS = [None]*count  
  for i in range (0, count):
    EXT_VARS[i] = data[i,:]  
#print EXT_VARS
t=EXT_VARS[1][0]
U=EXT_VARS[1][1]
beta=EXT_VARS[1][2]
mu=EXT_VARS[1][3]
nu_ext=(2*EXT_VARS[1][4]+1)*pi/beta
px_ext=EXT_VARS[1][5]
py_ext=EXT_VARS[1][6]
lat_const=EXT_VARS[1][7]
print
print "Evaluation for the following parameters:"
print 't = ', t
print 'U = ', U
print 'beta = ', beta
print 'mu = ', mu
print 'nu_ext = ', nu_ext
print 'px_ext = ', px_ext
print 'py_ext = ', py_ext
print 'a = ', lat_const
print
#End
###################################


p_old=[-1.5952651350881368, 2.1400079079679992]
p_new=[0.46625068823159177, -2.8691297014381987, -0.977206434208334, -0.08287685719422111]
print 'integrand_old = ', Integrand_1st_order(p_old)
print 'integrand_new = ', Integrand_2nd_order(p_new)

options={'epsrel':1.0e-6,'epsabs':1.0e-6,'limit':200}
myopts=[ options, options, options, options]


self_2nd_real,err_self_2nd_real=scipy.integrate.nquad(lambda p_x1,p_y1,p_x2,p_y2: scipy.real(integrand_real([p_x1,p_y1,p_x2,p_y2])),[[-pi,pi],[-pi,pi],[-pi,pi],[-pi,pi]], opts=myopts)
print '2nd: Real Part Det = ', self_2nd_real,err_self_2nd_real

self_2nd_imag,err_self_2nd_imag=scipy.integrate.nquad(lambda p_x1,p_y1,p_x2,p_y2: scipy.real(integrand_imag([p_x1,p_y1,p_x2,p_y2])),[[-pi,pi],[-pi,pi],[-pi,pi],[-pi,pi]], opts=myopts)
print '2nd: Imaginary Part Det = ', self_2nd_imag,err_self_2nd_imag

a2 = self_2nd_real

self_2nd_real,err_self_2nd_real=scipy.integrate.nquad(lambda p_x1,p_y1: scipy.real(integrand_real([p_x1,p_y1])),[[-pi,pi],[-pi,pi]], opts=myopts)
print 'First Order = ', self_2nd_real,err_self_2nd_real

print self_2nd_real + a2
