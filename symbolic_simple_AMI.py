'''
Symbolic AMI with simple poles
Updates: 2019-Feb-28
	 2019-March-1
         2019-March-2
         2019-March-3
         2019-March-4 
         2019-March-8
'''

from __future__ import division
from math import *
import math
import numpy
import numpy as np
import pylab as plt
import time
from copy import deepcopy
import copy
start_time = time.time()

##################################
###Definition of Fermi Function###
##################################
def fermi(e, m, B):  # Definition of fermi function
  val = 1.0/(cmath.exp(B*(e-m))+1.0)
  return val

#################################################
##Dispersion relation of 2D tight-binding model##
#################################################
def E( px, py, t, a): # 2D tight-binding dispersion
  val = 2*t*(cos(a*px) + cos(a*py)) 
  return val

'''
#################################
#####Symbolic Representation#####
#################################
def sym_rep_g(g):
  num_l=len(g)
  num_f=len(g[0])
  e=[]
  for j in range (0,num_l):
    a=[]
    a.append(j+1)
    a.append(1)
    for i in range (0,num_f):
      a.append(g[j][i])
    e.append(a)
  val=[]
  for i in range(0,num_l):
    val.append(e[i])
    val.append(g[i])
  
  return g,e,val
'''
#################################
#####Symbolic Representation#####
#################################
def sym_rep_g(g):  # symbolic representation of multiplication of Green's functions
  num_l=len(g)
  e=[[0]*(num_l)for i in range (0,num_l)]
  for j in range (0,num_l):
    e[j][j]=1
  val=[]
  for i in range(0,num_l):
    val.append(e[i])
    val.append(g[i])
  #print 'g=',g
  #print 'e=',e
  #print 'val=',val
  return g,e,val  # g:frequency part, e: momentum part, val: total representation 

################################
###########one sum##############
################################
# To do one summation over p-th frequency
def one_sum_n_G(G_fr, G_mo, p):  # G_fr: frequency part, G_mo: momentum part, p: frequency 0,1,...,n-1
  z = []
  g_freq = []  # To store freequency part
  g_mnta= []   # To store momentum part
  alpha_p = []
  #G_full = deepcopy(G_Init[2])
  G_freq = deepcopy(G_fr)  # Frequency part
  G_mnta = deepcopy(G_mo)  # Momentum part
  #g_full = deepcopy(G_Init[2])
  ###g_freq = deepcopy(G_Init[0])
  ###g_mnta = deepcopy(G_Init[1])
  
  val_freq = [[1]*(len(G_freq[0])) for i in range(len(G_freq))]  # To store freequency part
  val_mnta = [[1]*(len(G_mnta[0])) for i in range(len(G_mnta))]  # To store momentum part
  #for i in range (0, len(val_freq)):
    #for j in range (0,len(val_freq[i])):
      #val_freq[i][j]=1  
  #for i in range (0, len(val_mnta)):
    #for j in range (0,len(val_mnta[i])):
      #val_mnta[i][j]=1
  
  root_freq = []  # To store poles' frequency
  #root_mnta = []  # To store poles' momentum
  new_z=[]
  sign = []  # To store signs
  #root = []  # To store roots
  for i in range(0, len(G_freq)):
    g_freq.append(0)   # To store frequency part
    z.append(None)
    alpha_p.append(0)  # Initialize alpha values 
    root_freq.append(None)
    sign.append(1)
  for i in range(0, len(G_mnta)):
    g_mnta.append(0)
    #root_mnta.append(None)
  for i in range(0, len(G_freq)):  
      g_freq[i] = deepcopy(G_freq[i]) # Initialize g_freq
  for i in range(0, len(G_mnta)):     
      g_mnta[i] = deepcopy(G_mnta[i]) # Initialize g_mnta
  num_mnta=len(g_mnta)  # Dimension of g_mnta
  root_mnta=[[None]*(len(G_mnta[0]))for i in range (0,num_mnta)]  # Initialize root_mnta
  #print 'val_freq = ', val_freq
  #print 'val_mnta=', val_mnta
  #print 'g_full =', g_full
  #print 'g_freq =', g_freq
  #print 'g_mnta =', g_mnta
  #print 'root_mnta = ', root_mnta
  for i in range(0, len(g_freq)):  # Sweep Green's functions
    alpha_p[i] = g_freq[i][p]  # Compute alpha values
  #print 'alpha_p =', alpha_p
  index = -1
  for j in range (0, len(G_freq)):
    for i in range(0, len(G_freq)):  
      g_freq[i] = deepcopy(G_freq[i])  # Initialize g_freq
    for i in range(0, len(G_mnta)):    
      g_mnta[i] = deepcopy(G_mnta[i])  # Initialize g_mnta
    for i in range(0, len(G_freq)):
      alpha_p[i] = g_freq[i][p]  # Calculate alpha for pth frequency
    #if alpha_p[j]!=0:
      #for t in range (0,len(g_mnta)):
        #if g_mnta[t][0]==j+1:
          #g_mnta[t][1]=-1*alpha_p[j]*G_mnta[t][1]
          #root_mnta.append(g_mnta[t])
      #g_freq[j][p] = 0
      #root_freq.append(-alpha_p[j]*g_freq[j])  
    #print 'k=',k
    #print 'j=',j
    #print 'root_mnta = ', root_mnta  
    #print 'root_freq=', root_freq
    for k in range (0, len(G_freq)): # Sweep Green's functions 
        
           
        if k!=j:    # To make sure this is a different Green's functions      
          index = k   
          g_freq[index][p] = 0  # Set pth coefficient equal zero for others
          g_freq[j][p] = 0      # Set pth coefficient equal zero for jth Green's function
          #g_mnta[index][p] = 0
          #g_mnta[j][p] = 0
          #print g_freq
          if alpha_p[j]!=0:  # If pth frequency is present
            z[j] = -alpha_p[j]*g_freq[j][0]  # Get the zeroth element as flag
            root_freq[j]=-alpha_p[j]*g_freq[j] # Find frequency of the pole
            for i in range(0,len(g_mnta[j])): 
              root_mnta[j][i]=-alpha_p[j]*g_mnta[j][i]  # Find the momnetum of the pole
            
          for i in range(0, len(g_freq[j])):
            g_freq[index][i] = g_freq[index][i]-alpha_p[index]*alpha_p[j]*g_freq[j][i] # Find the residue: frequency part
          for i in range(0, len(g_mnta[j])):
            g_mnta[index][i] = g_mnta[index][i]-alpha_p[index]*alpha_p[j]*g_mnta[j][i] # Find the residue: momentum part
          #for i in range(0, len(g_mnta)):
            
          #if z[j] == 0:
           # z[j] = 'None'
          if alpha_p[j] == +1:
            sign[j] = +1  # Get sign
          if alpha_p[j] == -1:
            sign[j] = -1  # Get sign
          #else: 
          val_freq[j][index] = g_freq[index]  # Store residue: frequency part
          val_mnta[j][index] = g_mnta[index]  # Store residue: momentum part
          #print val[j][index], str(sign[j]), z[j]
          
          for i in range(0, len(G_freq)):  
            g_freq[i] = deepcopy(G_freq[i]) # Initialize g_freq
          for i in range(0, len(G_mnta)):    
            g_mnta[i] = deepcopy(G_mnta[i])  # Initialize g_mnta  
          for i in range(0, len(G_freq)):
            alpha_p[i] = g_freq[i][p]       
  #print 'root_freq=', root_freq
  #print 'root_mnta=', root_mnta
  #print 'val_freq=', val_freq
  #print 'val_mnta=', val_mnta
  #print G_freq
  #print G_mnta
  #print 'z=',z
  roots_number = 0  # To find the number of simple poles
  for i in range (0, len(z)):
    if z[i] != None:  # If there is a pole
      roots_number = roots_number + 1
  j = 0
  new_val_freq = []  # To store residue: frequency part
  new_val_mnta = []  # To store residue: momentum part
  new_sign = []      # To store sugn
  new_root_freq = [] # To store poles: frequency part
  new_root_mnta = [] # To store poles: momentum part
  new_z=[]
  #fin_val = []
  for i in range(0, roots_number):
    new_z.append(None)
    new_root_freq.append(0)
    new_root_mnta.append(0)
    new_val_freq.append(0)
    new_val_mnta.append(0)
    new_sign.append(0)
    #fin_val.append(0)
  #print "roots_number = ", roots_number
  for i in range (0, len(z)): 
    if z[i] != None:
      new_root_freq[j] = root_freq[i]
      new_root_mnta[j] = root_mnta[i]
      new_z[j] = z[i]
      new_val_freq[j] = val_freq[i]
      new_val_mnta[j] = val_mnta[i]
      new_sign[j] = sign[i]
      j=j+1  
  #print 'new_root_freq = ', new_root_freq
  #print 'new_root_mnta = ', new_root_mnta
  #print 'new_z = ', new_z
  #print 'new_val_freq = ', new_val_freq
  #print 'new_val_mnta = ', new_val_mnta
  #print j
  #print val[0][0]
  for i in range (0, len(G_freq)):
    k = 0
    for l in range (0, len(G_freq)):
      if type(val_freq[i][l]).__name__ != 'int':
        k=k+1 
  #print k
  fin_val_freq = [[1]*(k) for i in range(0, roots_number)]
  fin_val_mnta = [[1]*(k) for i in range(0, roots_number)]
  #print fin_val
  #print new_val[0][1]
  s = 0
  for i in range (0, roots_number):
    s = 0
    for l in range (0, len(G_freq)):
      if type(new_val_freq[i][l]).__name__ != 'int': 
        #print 'i = ', i, 'l = ', l, new_val[i][l]
        fin_val_freq[i][s] = new_val_freq[i][l]
        fin_val_mnta[i][s] = new_val_mnta[i][l]
        s = s + 1
  #print 'new_sign = ', new_sign
  #print 'new_root_freq = ', new_root_freq
  #print 'new_root_mnta = ', new_root_mnta  
  #print 'fin_val_freq = ', fin_val_freq
  #print 'fin_val_mnta = ', fin_val_mnta
  #print 'new_z = ', new_z 
  #print 'Done! See G_out.txt' 
  return new_sign, new_root_freq, new_root_mnta, fin_val_freq, fin_val_mnta, new_z
###################################

###################################
#########Initial energies##########
###################################
# To evaluate the numerical values of initial dispersions
def Energy(g_freq_init,p,h,M):
  l = len(g_freq_init)
  e=[None] *l
  
  for i in range(0,l):
    val=[0]*2
    for j in range (0,len(g_freq_init[0])):
      
      for r in range(0,len(p[0])):
        val[r] = p[j][r]*g_freq_init[i][j]+val[r]
    
    e[i]=-2*h*(cos(val[0])+cos(val[1]))-M 
   
  return e

###################################
#########Energy evaluation#########
###################################
# To evaluate the numerical values of dispersions
def Energy_eval(e_array,g_freq_init,p,h,M):
  e=Energy(g_freq_init,p,h,M)
  val=0
  for i in range(0,len(e_array)):
    val = val - e_array[i]*e[i]
  #print e  
  return val

###################################
###########Poles merging###########
###################################
#To merge frequency and momenta poles
def Poles(p_f, p_m, g_freq_init,p,h,M): # p_f: frequency array, p_m: momenta array of poles
 
  num_len = len(g_freq_init[0])+1
  p_f_copy = deepcopy(p_f)
  term = [None]*num_len
  for i in range (0,len(p_f)):
    for j in range (0,len(p_f[i])):
      if p_f[i][j]!=None:
        #print 'p_f = ', p_f[i][j]
        #print 'p_m = ', p_m[i][j]
        for k in range (0,len(p_f[i][j])):
          #print 'p_f[i][j][k] = ', p_f[i][j][k]  
          #print 'p_m[i][j][k] = ', p_m[i][j][k]
          term[0] = Energy_eval(p_m[i][j][k],g_freq_init,p,h,M)
          #print 'term[0] = ', term[0] 
          for r in range (1, num_len):
            term[r] = p_f[i][j][k][r-1]
          #print 'term = ', term
          p_f[i][j][k] = deepcopy(term)
  return p_f


###################################
#########fermi/bose function#######
###################################
# Definition of f Function
def f(z,B):
  lz = len(z)
  index = 0
  nu = [0]*lz
  for i in range (1, lz):
    if (z[i]!=0):
      index = index + 1
  sigma = (-1)**index
  #print index
  val = 1.0/(sigma*exp(B*z[0]-0.01)+1.0)
  return val
###################################  

###################################
############dot_num################
###################################
# Definition of dot operation between two array of numbers
def dot_num(a, b):
  la = len(a)
  lb = len(b)
  if la!=lb:
    print 'BEEP! Two arrays mush have the same length!'
  else:
    l = la
    val = [None]*l
    for i in range(0, l):
      val[i] = a[i]*b[i]
    #print 'A dot B = ', a, '.', b, '= ', val
    #return 'Done!' 
  return val

###################################
############dot_arr################
###################################
# Definition of dot operation between array of numbers and array of arrays of numbers
def dot_arr(c, d):
  lc = len(c)
  ld = len(d)
  count = 0
  if lc!=ld:
    print 'BEEP! Two arrays must have the same length!'
  if lc==ld:
    l = lc
    result = [None]*l
    for i in range(0,l):
      
      if len(c[i]) == len (d[i]):
        result[count] = dot_num(c[i],d[i])
        count = count + 1
      else:
        print 'BEEP! Two arrays elements must have the same length!'
        return None
  return result   

###################################
##############cross################
###################################
# Definition of cross operation  
def cross(a, c):
  la = len(a)
  lc = len(c)
  lval = 0
  count = 0
  if la!=lc:
    print 'BEEP! Two arrays must have the same length!'   
  if la == lc:
    l = la
    for i in range (0, l):
      lci = len(c[i])
      for j in range (0, lci):
        lval = lval + 1
    val = [None]*lval
    for i in range (0, l):
      lci = len(c[i])
      for j in range (0, lci):
        val[count] = a[i] * c[i][j]
        count = count + 1
    #print count
    #print lval
    #print val
    return val 

###################################
############Evaluation#############
###################################
def G_eval(g_arr_freq, g_arr_mnta, ext, g_freq_init,p,h,M,bet):
  val = 1
  delta = 0.1
  num = len(g_arr_freq[0])
  #print 'num = ',num
  for i in range (0, len(g_arr_freq)):
    e=Energy_eval(g_arr_mnta[i],g_freq_init,p,h,M)
    #print e
    if (e + g_arr_freq[i][num-1]*1j*ext )!=0:
      val = val/( e + g_arr_freq[i][num-1]*(ext*1j) )
    else:
      val = 0
  return val






