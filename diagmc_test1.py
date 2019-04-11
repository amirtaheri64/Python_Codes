'''
Diagrammatic Monte Carlo test for

1 + \int_a^b dxexp(-x) + \int_a^b dy \int_a^b dz exp(-y-z)

'''

###################################
##########Import libraries#########
###################################
#Begin
from math import *
from random import random
from random import randint
import numpy

#End
###################################

###################################
#########First integrand###########
###################################
#Begin
def first_integrand(X):
  return exp(-X)
#End
###################################

###################################
#########Second integrand##########
###################################
#Begin
def second_integrand(Y,Z):
  return exp(-Y-Z)
#End
###################################

###################################
############Integrand##############
###################################
#Begin
def integrand(ord_var):  # Return D functions
  l=len(ord_var)
  s=0
  for i in range (1,l):
    s=s+ord_var[i]  
  return exp(-s)
#End
###################################

###################################
######Proposal probabilities#######
###################################
#Begin
def W(ORDER_NEW,ORDER_OLD): 
  val=1.0/(b-a)**ORDER_NEW
  #if ORDER_NEW==0:
    #val = 1.0
  #if ORDER_NEW==1:
    #val = 1.0/(b-a)
  #if ORDER_NEW==2:
    #val = 1.0/(b-a)/(b-a)  
  return val
#End
###################################

###################################
###########Update_order############
###################################
#Begin
def update_order(ORDER):  # Change order by 1
  r=random()
  if r<1.0/2.0: 
    val=ORDER-1
  if r>=1.0/2.0:
    val=ORDER+1
  #val=randint(0,1)
  return val
#End
###################################

###################################
#####Update internal variables#####
###################################
#Begin
def update_int_var(ORDER):  # Update internal variables
  val=[]
  val.append(ORDER)
  if ORDER==0 or ORDER==-1 or ORDER==3:
    return val
  if ORDER==1:
    X=(b-a)*random()+a
    val.append(X)
    return val
  if ORDER==2:
    Y=(b-a)*random()+a
    val.append(Y)
    Z=(b-a)*random()+a
    val.append(Z)
    return val
#End
###################################
    

###################################
#########Acceptance ratio##########
###################################
#Begin
def acp_rt(ord_var_old,ord_var_new):  # Return the acceptance ratio
  #print 'old state = ', ord_var_old
  #print 'new state = ', ord_var_new
  ORDER_OLD=ord_var_old[0]
  ORDER_NEW=ord_var_new[0]
  if ORDER_NEW==3 or ORDER_NEW==-1:  # If the order is out of range
    val=0
    return val
  val=integrand(ord_var_new)*W(ORDER_OLD,ORDER_NEW)/integrand(ord_var_old)/W(ORDER_NEW,ORDER_OLD)
  return val
#End
###################################

###################################
########Collect statistics#########
###################################
#Begin
def new_state(OLD_STATE):  # Accept or reject the proposed state
  PROPOSED_ORDER = update_order(OLD_STATE[0])
  #print 'order_old = ', ORDER_OLD
  #print 'proposed order = ', PROPOSED_ORDER	
  PROPOSED_STATE = update_int_var(PROPOSED_ORDER)
  A=acp_rt(OLD_STATE,PROPOSED_STATE)
  #print 'A = ', A
  r=random()
  if A <= r:
    NEW_STATE = OLD_STATE
  else:  
    NEW_STATE = PROPOSED_STATE
  return NEW_STATE
#End
###################################


a=0  # Lower limit of the integral
b=2  # Upper limit of the integral
x=[]  # To store Markov's chain
var=0
avg_x=0

old_order=0  # Start from the zeroth order
old_state=update_int_var(old_order)[:]  # Generate the initial state
N=10000 # Length of Markov chain
MC_iter=100000 # Monte Carlo iteration
avg=0  # To store average
for j in range (0,MC_iter):
  N_0=0.0
  for i in range (0,N):  # To propose N states
    next_state = new_state(old_state)[:]  # Update 
    old_state = next_state[:]  
    if old_state[0] == 0: # To compute N_0
      N_0 = N_0 + 1
  x.append(N_0)
  avg=avg+N_0 # To compute the average of N_0
  #print N_0
N_0=avg/MC_iter  
print 'a = ', a
print 'b = ', b
print 'N = ', N
print 'MC_iter = ', MC_iter
for i in range (0,len(x)):
  avg_x = avg_x + x[i]
avg_x=avg_x/len(x)
print 'avg_x= ', avg_x
for i in range (0,MC_iter):  # Calculate the variance
  var = var + (x[i]-avg_x)**2
var = sqrt(var/MC_iter-1)
var = var/sqrt(MC_iter)
print 'uncertainy = ', N*var/avg_x/avg_x
print 1+(N-N_0)/N_0
var=numpy.std(x)
var = var/sqrt(MC_iter)
print 'uncertainy = ', N*var/avg_x/avg_x



