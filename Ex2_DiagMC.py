'''
DiagMC: Second Example: computation of 1+int_a^b xdx
Date 15-NOV-2018
'''

##################################
##########Import Libraries########
##################################
#Begin
from random import *
from math import *
#End
###################################

###################################
######Proposal probabilities#######
###################################
#Begin
def W(n2,n1,a,b):  # n: D function order, a and b: Integral limits
  if n2==0:
    val = 1.0/2.0
  if n2==1:
    val = 1.0/2.0/(b-a)
  return val
#End
###################################

###################################
####Definition of D function#####
###################################
#Begin
def D(n,x): # n: D Function order, x: D function argument
  if n==0:
    return 1.0
  else:
    return x
#End
###################################


###################################
##########Minimum finder###########
###################################
#Begin
def min(a,b):
  if a<=b:
    return a
  else:
    return b
#End
###################################

###################################
#######Accept-Reject ratios########
###################################
#Begin
def A(n2,n1,x2,x1,a,b):  # n: D Function order, x: D function argument, x: D function argument
  val = D(n2,x2)*W(n1,n2,a,b)/D(n1,x1)/W(n2,n1,a,b)
  return min(1.0,val)
#End
###################################

###################################
######Propose the next state#######
###################################
#Begin
def propose_state(a,b):
  r = random()
  if r<=0.5:  # Go to the zeroth order state
    return 0, 1   
  else:       # Go to the first order state
    r = random()
    return 1, r*(b-a)+a

#End
###################################

###################################
###########Update state############
###################################
#Begin
def update(CURRENT_STATE,a,b):  # a and b: Integral limits
  new_state= [None]*2
  PROPOSED_STATE = propose_state(a,b)
  r = random()
  #print 'proposed = ', PROPOSED_STATE
  #print 'r = ', r
  #print 'A = ', A(PROPOSED_STATE[0],CURRENT_STATE[0],PROPOSED_STATE[1],CURRENT_STATE[1],a,b)
  if A(PROPOSED_STATE[0],CURRENT_STATE[0],PROPOSED_STATE[1],CURRENT_STATE[1],a,b) <= r:
    new_state = CURRENT_STATE
  else:  
    new_state = PROPOSED_STATE
  return new_state
#End
###################################

A_l = 2 # Define a
B_l = 5 # Define b
N = 10000 # Total sequence length
MC_Iter = 100000  # Number of Iterations
N_0 = 0.0
avg = 0.0
x=[]
y=[]
var = 0

###################################
###########Check Point#############
###################################
#Begin 
#print 'W = ' , W(1,1,A_l,B_l)   # OK
#print 'D = ', D(1,2.0)  # OK
#print 'A(x",x) = ', A(1,1,5,4,A_l,B_l) # OK
#current_state = propose_state(A_l,B_l)  # OK
#print current_state  # OK
#print update(current_state,A_l,B_l)  # OK
#End
###################################


current_state = propose_state(A_l,B_l)
for j in range (0,MC_Iter):
  for i in range (0,N):
    next_state = update(current_state,A_l,B_l)
    current_state = next_state
    if current_state[0] == 0:
      N_0 = N_0 + 1
  avg = avg + N_0
  x.append(N_0)
  #print N_0
  #print x[j]
  N_0 = 0
N_0 = avg/MC_Iter
print float(N)/N_0
avg_x=0
for i in range (0,len(x)):
  avg_x = avg_x + x[i]
avg_x=avg_x/len(x)
print 'avg_x= ', avg_x
print float(N)/avg_x 


for i in range (0,MC_Iter):  # Calculate the variance
  var = var + (x[i]-avg_x)**2
var = sqrt(var/MC_Iter)
var = var/sqrt(MC_Iter)
print 'uncertainy = ', N*var/avg_x/avg_x
#print 'y = ', y 
'''
for i in range (0,len(y)):  # Calculate the variance
  var = var + (y[i]-avg_y)**2
var = sqrt(var/len(y))
print 'std_err = ', var/sqrt(len(y)) 
'''
'''
print 'N = ', N, 'T = ', MC_Iter, 'N_0 = ', N_0

avg = 0
for i in range (0,MC_Iter):  # Calculate the average
  avg = avg + x[i]
#avg = avg/MC_Iter
print 'avg = ', avg
for i in range (0,MC_Iter):  # Calculate the variance
  var = var + (x[i]-avg)**2
var = sqrt(var/MC_Iter)

print 'Result = ', 1 + (N-N_0)/N_0
print 1+(N-N_0)/N_0
print 'std_err = ', var/sqrt(MC_Iter)  
'''

