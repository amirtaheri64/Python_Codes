'''
DiagMC: First Example: computation of 1+C
Date 15-NOV-2018
'''

##################################
##########Import Libraries########
##################################
#Begin
from random import *
#End
###################################

###################################
######Proposal probabilities#######
###################################
#Begin
def W(x2,x1,c):
  if (x1==c and x2==c) or (x1==1 and x2==1) or (x1==1 and x2==c) or (x1==c and x2==1):
    val = 1.0/2.0
    return val
  else:
    print 'OOPS!'
#End
###################################

###################################
#####Definition of D function######
###################################
#Begin
def D(x,C):
  if x==1 or x==C:
    return x 
  else:
    print 'OOPS!'
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
def A(x2,x1,c):
  val = D(x2,c)*W(x1,x2,c)/D(x1,c)/W(x2,x1,c)
  return min(1.0,val)
#End
###################################

###################################
######Propose the next state#######
###################################
#Begin
def propose_state(c):
  r = random()
  if r<=0.5:
    return 1
  else:
    return c
#End
###################################

###################################
###########Update state############
###################################
#Begin
def update(CURRENT_STATE,c):
  PROPOSED_STATE = propose_state(c)
  a = random()
  if A(PROPOSED_STATE,CURRENT_STATE,c) <= a:
    new_state = CURRENT_STATE
  else:  
    new_state = PROPOSED_STATE
  return new_state
#End
###################################

C=4  # Define C
N = 5000  # Total sequence length
MC_Iter = 100  # Number of Iterations
N_0 = 0.0
avg = 0.0

###################################
###########Check Point#############
###################################
#Begin 
#print 'W = ' , W(1,1,C)
#print 'D = ', D(2,C) 
#print 'A(1,1) = ', A(1,1,C)
#print 'A(1,C) = ', A(1,C,C)
#print 'A(C,1) = ', A(C,1,C)
#print 'A(C,C) = ', A(C,C,C)
#End
###################################

current_state = propose_state(C)
for j in range (0,MC_Iter):
  for i in range (0,N):
    next_state = update(current_state,C)
    current_state = next_state
    if current_state == 1:
      N_0 = N_0 + 1
  avg = avg + N_0
  N_0 = 0
N_0 = avg/MC_Iter
print 'N = ', N, 'N_0 = ', N_0
print '1+C = ', 1 + (N-N_0)/N_0

  


