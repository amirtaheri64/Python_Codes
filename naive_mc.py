'''
Monte Carlo Integration: naive computation

#########################################################
Check points: 

From explicit integrands:

Hartree =  (-0.83626183323827918, 0.00035564681976439974)

Second order term: Imag part =  -0.0898533208632 0.0002745548805
Second order term: Real part =  0.221504891341 0.0011175190974


#########################################################
'''

###################################
########Import Libraries###########
###################################
#Begin
from __future__ import division
from skmonaco import mcquad
from math import *
import scipy.integrate 
import numpy
import time
import numpy as np
import pylab as pl
from random import random
start_time = time.time()
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
  val = -2*t*(cos(p[0]*a)+cos(p[1]*a))-mu
  return val
#End
###################################

###################################
########Integrand Hartree##########
###################################
#Begin
def Integrand_1st_order_new(p): # p: momentum, m: chemical potential, t: hopping amplitude, a: lattice constant, B: inverse temperature, u: Hubbard potential
  dis_en=E(p)
  val=-U*fermi(dis_en)/2/2/pi/pi
  return val
###################################

###################################
########Pick a random point########
###################################
#Begin
def pick(a,b):
  r1=random()
  new=r1*(b-a)+a
  return new 
#End
###################################

###################################
#########Naive MC Hartree##########
###################################
#Begin
def naive_mc_Hartree(dim,N):
  hits=0
  p=[None]*dim
  chain=[]
  for i in range (0,N):
    for j in range (0,dim):
      p[j]=pick(-pi,pi)  # Pick momentum
      
    hits=hits+Integrand_1st_order_new(p)
    chain.append(Integrand_1st_order_new(p))
  
  s=0
  avg=hits/N  
  measurement=((2*pi)**dim)*hits/N   
  for i in range (0,N):
    s=s+(chain[i]-avg)*(chain[i]-avg)
  s=s/N
  s=sqrt(s)
  uncertainty = s/sqrt(N)*(2*pi)**dim 
  return measurement, uncertainty
#End
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

############Integrand##############
###################################
#Begin
def Integrand_3rd_order(p): # p: momentum, m: chemical potential, t: hopping amplitude, a: lattice constant, B: inverse temperature, u: Hubbard potential
  k=[p[0],p[1]]
  e1=E(k)
  k = [p[2],p[3]]
  e3=E(k)
  k = [p[4],p[5]]
  e4=E(k)
  k=[px_ext+p[2]-p[0],py_ext+p[3]-p[1]]
  e5=E(k)
  val1=fermi(e1)*fermi(e3)*fermi(e4)
  val2=fermi(e1)*bose(e1+e5)*fermi(e4)
  val3=-fermi(-e5)*bose(e1+e5)*fermi(e4)
  val4=-fermi(-e5)*fermi(e3)*fermi(e4)
  val=val1+val2+val3+val4
  val=-U*U*U*val/2/2/2/2/2/2/pi/pi/pi/pi/pi/pi  # Only the prefactor
  val=val/(1j*nu_ext-e1-e5+e3)/(1j*nu_ext-e1-e5+e3)
  return val
#End
###################################

###################################
#######Integrand real part#########
###################################
#Begin
def Int_2nd_Real(p):
  val = float( Integrand_2nd_order(p).real )
  return val
#End
###################################

###################################
#######Integrand imag part#########
###################################
#Begin
def Int_2nd_Imag(p):
  val = float( Integrand_2nd_order(p).imag )
  return val
#End
###################################

###################################
#######Integrand real part#########
###################################
#Begin
def Int_3rd_Real(p):
  val = float( Integrand_3rd_order(p).real )
  return val
#End
###################################

###################################
#######Integrand imag part#########
###################################
#Begin
def Int_3rd_Imag(p):
  val = float( Integrand_3rd_order(p).imag )
  return val
#End
###################################


###################################
#############Naive MC##############
###################################
#Begin
def naive_mc_2nd(dim,N):
  hits_imag=0
  hits_real=0
  p=[None]*dim
  chain_imag=[]
  chain_real=[]
  for i in range (0,N):
    for j in range (0,dim):
      p[j]=pick(-pi,pi)  # Pick momentum
      
    hits_imag=hits_imag+Int_2nd_Imag(p)
    chain_imag.append(Int_2nd_Imag(p))
    hits_real=hits_real+Int_2nd_Real(p)
    chain_real.append(Int_2nd_Real(p))
  
  s_imag=0
  avg_imag=hits_imag/N  
  #print 'avg =', avg
  #print chain
  measurement_imag=((2*pi)**dim)*hits_imag/N   
  for i in range (0,N):
    s_imag=s_imag+(chain_imag[i]-avg_imag)*(chain_imag[i]-avg_imag)
  s_imag=s_imag/N
  s_imag=sqrt(s_imag)
  uncertainty_imag = s_imag/sqrt(N)*(2*pi)**dim 

  s_real=0
  avg_real=hits_real/N  
  #print 'avg =', avg
  #print chain
  measurement_real=((2*pi)**dim)*hits_real/N   
  for i in range (0,N):
    s_real=s_real+(chain_real[i]-avg_real)*(chain_real[i]-avg_real)
  s_real=s_real/N
  s_real=sqrt(s_real)
  uncertainty_real = s_real/sqrt(N)*(2*pi)**dim 
  return measurement_imag, uncertainty_imag, measurement_real, uncertainty_real  
#End
###################################


###################################
######Read external variables######
###################################
'''
data = np.loadtxt('ext_vars.dat')  # reading external variables from a file
count = len(open('ext_vars.dat').readlines())
if count==1:
  EXT_VARS = [None]*count
  EXT_VARS[0]=data
if count>1:
  EXT_VARS = [None]*count  
  for i in range (0, count):
    EXT_VARS[i] = data[i,:]  
'''

#print EXT_VARS
#t=EXT_VARS[1][0]
#U=EXT_VARS[1][1]
#beta=EXT_VARS[1][2]
#mu=EXT_VARS[1][3]
#nu_ext=(2*EXT_VARS[1][4]+1)*pi/beta
#px_ext=EXT_VARS[1][5]
#py_ext=EXT_VARS[1][6]
#a=EXT_VARS[1][7]



t=1
a=1
U=1
beta=5
mu=0
nu_ext=pi/beta
px_ext=pi
py_ext=pi/3.0
kx1=ky1=0
kx2=pi
ky2=pi/3.0
'''
print
print "Evaluation for the following parameters:"
print 't = ', t
print 'U = ', U
print 'beta = ', beta
print 'mu = ', mu
print 'nu_ext = ', nu_ext
print 'px_ext = ', px_ext
print 'py_ext = ', py_ext
print 'a = ', a
print
k_init=[ kx1,ky1,kx2,ky2 ]
print k_init
print Integrand_2nd_order(k_init)

'''
'''
m=1    # order of diagram
DIM=m*2  # number of independent momentum variables: px and py
step=1000000
print 'Hartree = ', naive_mc_Hartree(DIM,step)
print
'''

step=10000000
m=2    # order of diagram
DIM=m*2  # number of independent momentum variables: px and py

print
print "Evaluation for the following parameters:"
print 't = ', t
print 'U = ', U
print 'beta = ', beta
print 'mu = ', mu
print 'nu_ext = ', nu_ext
print 'px_ext = ', px_ext
print 'py_ext = ', py_ext
print 'a = ', a
print
for i in range(1,10):
  print 'samples = ', step*(i+1)
  out_second_order = naive_mc_2nd(DIM,step*(i+1))
  
  print 'Second order term: Imag part = ', out_second_order[0], out_second_order[1]
  print 'Second order term: Real part = ', out_second_order[2], out_second_order[3]
  print

'''
k_init=[0]*6
print k_init
print Integrand_2nd_order(k_init)
'''
'''
m=1    # order of diagram
DIM=m*2  # number of independent momentum variables: px and py
step=1000000
print 'Hartree = ', naive_mc_Hartree(DIM,step)
print
'''
'''
inc=10000000
m=3    # order of diagram
DIM=m*2  # number of independent momentum variables: px and py


  

print
print "Evaluation for the following parameters:"
print 't = ', t
print 'U = ', U
print 'beta = ', beta
print 'mu = ', mu
print 'nu_ext = ', nu_ext
print 'px_ext = ', px_ext
print 'py_ext = ', py_ext
print 'a = ', a
print
for i in range(0,10):
  step = inc*(i+1)
  out_third_order = naive_mc_2nd(DIM,step)
  print 'Second order term: Imag part = ', out_third_order[0], out_third_order[1]
  print 'Second order term: Real part = ', out_third_order[2], out_third_order[3]
  #print
'''
