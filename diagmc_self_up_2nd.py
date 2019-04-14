'''
Diagrammatic Monte Carlo test for
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
######Proposal probabilities#######
###################################
#Begin
def W(ORDER_NEW,ORDER_OLD): 
  val=1.0/(b-a)**(2*ORDER_NEW)
  return val
#End
###################################

###################################
###########Update_order############
###################################
#Begin
def update_order(ORDER):  # Change order by 1
  r=random()
  if flag_re==1:
    if r<1.0/2.0: 
      val=ORDER-1
    if r>=1.0/2.0:
      val=ORDER+1
    return val
  if flag_re==0:
    if r<1.0/2.0: 
      val=ORDER-2
    if r>=1.0/2.0:
      val=ORDER+2
    return val
#End
###################################

###################################
#####Update internal variables#####
###################################
#Begin
def update_int_var(ORDER):  # Update internal variables
  val=[]
  for i in range (0,2*ORDER):
    X=(b-a)*random()+a
    val.append(X)
  return val
#End
###################################

###################################
#########Acceptance ratio##########
###################################
#Begin
def acp_rt(P_OLD,P_NEW):  # Return the acceptance ratio
  #print 'old state = ', ord_var_old
  #print 'new state = ', ord_var_new
  ORDER_OLD=len(P_OLD)/2
  ORDER_NEW=len(P_NEW)/2
  if ORDER_NEW>=3 or ORDER_NEW<=-1:  # If the order is out of range
    val=0
    return val
  if flag_re==1:
    if integrand_real(P_OLD)==0:
      val=1
      return val
    val=abs(integrand_real(P_NEW))*W(ORDER_OLD,ORDER_NEW)/abs(integrand_real(P_OLD))/W(ORDER_NEW,ORDER_OLD)
  if flag_re==0:
    if integrand_imag(P_OLD)==0:
      val=1
      return val
    val=abs(integrand_imag(P_NEW))*W(ORDER_OLD,ORDER_NEW)/abs(integrand_imag(P_OLD))/W(ORDER_NEW,ORDER_OLD)
  #print 'integrand_new = ', integrand_real(P_NEW)
  #print 'integrand_old = ', integrand_real(P_OLD)
  return val
#End
###################################

###################################
###########Updatefunction##########
###################################
#Begin
def new_state(OLD_STATE):  # Accept or reject the proposed state
  ORDER_OLD=len(OLD_STATE)/2
  PROPOSED_ORDER = update_order(ORDER_OLD)
  #print 'order_old = ', ORDER_OLD
  #print 'proposed order = ', PROPOSED_ORDER	
  PROPOSED_STATE = update_int_var(PROPOSED_ORDER)
  #print 'proposed state = ', PROPOSED_STATE
  
  R=acp_rt(OLD_STATE,PROPOSED_STATE)
  #print 'R = ', R
  r=random()
  #print 'r = ', r
  if R <= r:
    NEW_STATE = OLD_STATE
  else:  
    NEW_STATE = PROPOSED_STATE
  if flag_re==1:
    if integrand_real(NEW_STATE)>0:
      sign=1
    else:
      sign=-1
  if flag_re==0:
    if integrand_imag(NEW_STATE)>0:
      sign=1
    else:
      sign=-1
  #print 'state_new = ', NEW_STATE
  return NEW_STATE,sign
#End
###################################

###################################
############Measurement############
###################################
#Begin
def measure(OLD_STATE):
  x0=[]  # To store Markov's chain
  x1=[]  # To store Markov's chain
  x2=[]  # To store Markov's chain
  for j in range (0,MC_iter):
    N0=0.0
    N1=0.0
    N2=0.0
    for i in range (0,N):  # To propose N states
      next_state = new_state(OLD_STATE)[:]  # Update 
      OLD_STATE = next_state[:][0]  
      if i>1000:
        if len(OLD_STATE)/2 == 0: # To compute N0
          N0 = N0 + next_state[1]
        if len(OLD_STATE)/2 == 1: # To compute N1
          N1 = N1 + next_state[1]
        if len(OLD_STATE)/2 == 2: # To compute N2
          N2 = N2 + next_state[1]
    x0.append(N0)
    x1.append(N1)
    x2.append(N2)
  return x0,x1,x2
#End
###################################

###################################
###############Error###############
###################################
#Begin
def stats(x0,x1,x2):  # xi: Markov's chain corresponding to Ni
  var0=0.0  
  var1=0.0
  var2=0.0
  avg_0=0.0
  avg_1=0.0
  avg_2=0.0
  error_N1_N0=0.0
  error_N2_N0=0.0

  for i in range (0,len(x0)):
    avg_0 = avg_0 + x0[i]
  avg_0=avg_0/len(x0)
  for i in range (0,len(x1)):
    avg_1 = avg_1 + x1[i]
  avg_1=avg_1/len(x1)
  for i in range (0,len(x2)):
    avg_2 = avg_2 + x2[i]
  avg_2=avg_2/len(x2)
  #print 'avg_0= ', avg_0
  #print 'avg_1= ', avg_1
  #print 'avg_2= ', avg_2

  #var0=np.std(x0)
  #var0 = var0/sqrt(MC_iter)
  #var1=np.std(x1)
  #var1 = var1/sqrt(MC_iter)
  #var2=np.std(x2)
  #var2 = var0/sqrt(MC_iter)
  for i in range (0,MC_iter):  # Calculate the variance
    var0 = var0 + (x0[i]-avg_0)**2
    var1 = var1 + (x1[i]-avg_1)**2
    var2 = var2 + (x2[i]-avg_2)**2
  var0 = sqrt(var0/(MC_iter-1))
  var0 = var0/sqrt(MC_iter)
  #print 'error_0 = ', var0
  var1 = sqrt(var1/(MC_iter-1))
  var1 = var1/sqrt(MC_iter)
  #print 'error_1 = ', var1
  var2 = sqrt(var2/(MC_iter-1))
  var2 = var2/sqrt(MC_iter)
  #print 'error_2 = ', var2

  a1 = avg_1/avg_0
  #print 'N1/N0 = ', a1
  if a1!=0:
    #error_N1_N0 = abs(a1) * sqrt( (var0/avg_0)**2 + (var1/avg_1)**2 )
    error_N1_N0 = abs(a1) * (var0/abs(avg_0) + var1/abs(avg_1) )
    #print 'error_N1/N0 = ', error_N1_N0 
  a2=avg_2/avg_0
  #print 'N_2/N_0 = ', a2
  if a2!=0:
    #error_N2_N0 = abs(a2) * sqrt( (var0/avg_0)**2 + (var2/avg_2)**2 )
    error_N2_N0 = abs(a2) * ( var0/abs(avg_0) + var2/abs(avg_2) )
    #print 'error_N2/N0 = ', error_N2_N0
  Result = a1+a2
  #print 'Result = ', Result
  if avg_1!=0:
    #error = sqrt (error_N1_N0**2 + error_N2_N0**2)
    error = error_N1_N0 + error_N2_N0
  else:
    error = error_N2_N0
  #print 'error = ', error
  return a1, error_N1_N0, a2, error_N2_N0, a1+a2, error
#End
###################################


###################################
###############Main################
###################################
#Begin

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

flag_re = 0  # Global variable  
a=-pi  # Lower limit of the integrals
b=pi   # Upper limit of the integrals
old_order=0  # Start from the zeroth order
old_state=update_int_var(old_order)[:]  # Generate the initial state
N=1000000 # Length of Markov chain
MC_iter=1000 # Monte Carlo iteration
#print 'a = ', a
#print 'b = ', b
print 'N = ', N
print 'MC_iter = ', MC_iter

chains = measure(old_state)  # Measurement
output = stats(chains[0],chains[1],chains[2])  # Output

if flag_re==1:
  print 'Results for the real part:'
else:
  print 'Results for the imaginary part:'

print 'a1 = ', output[0], '+-', output[1]
print 'a2 = ', output[2], '+-', output[3]
print 'Result = ', output[4], '+-', output[5]


#var=numpy.std(x)
#var = var/sqrt(MC_iter)
#print 'uncertainy = ', N*var/avg_x/avg_x
