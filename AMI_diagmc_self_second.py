'''
AMI+diagMC for the self energy function
up to the second order.

Last updates: APRIL-17-2019
'''

###################################
#########Import libraries##########
###################################
#Begin
from AMI_zero import *
from add_remove import *
from Hubbard_like_diagram import *
from is_connected import *
from Label_self import *
from Symbolic_multi_AMI_new import *
#End
###################################

###################################
###########Diagram checks##########
###################################
#Begin
def checks(G,M):
  FLAG_AMI=0
  if is_connected(G) and Irreducible(G):
    labels = G.es["Label"]
    reset_g(G,M)
    if Hubbard_diagram(G,M):
      AMI_INPUT = AMI_Input(labels,M)
      FLAG_AMI=AMI_zero(AMI_INPUT)
  val=[]
  val.append(FLAG_AMI)
  if FLAG_AMI!=0:
    val.append(AMI_INPUT)
  return val
#End
###################################

###################################
############Update momenta#########
###################################
#Begin
def update_mnta(G,M):
  out=[]
  ###################################
  #####Initialize moemntum array#####
  ###################################
  #Begin
  k_init=[ [0]*2 for i in range (0,M+1) ]  # initialize momenta
  k_init[M][0] = px_ext
  k_init[M][1] = py_ext
  #print k_init
  #End
  ###################################

  #start_time=time.time()
  #print 'check G = ', G
  #print 'order = ', M
  CHECK_diag = checks(G,M)
  flag_AMI = CHECK_diag[0]
  #out.append(flag_AMI)
  if flag_AMI==0 or M>=3 or M<=0:
    out.append(0)
    out.append(0)
    return out
  if flag_AMI==1:
    if M==1:
      out.append(1.0)
      out.append(1.0)
      return out
    ami_input=CHECK_diag[1]
    b= AMI_arrays_out(ami_input,M)   # Generate S, P, and R arrays and store in b
    #print("--- Construction_time ---", (time.time() - start_time))  # Compute the construction time

    ###################################
    #########Store AMI arrays########## 
    ###################################
    #Begin
    S_list=deepcopy(b[0])
    P_list_freq=deepcopy(b[1])
    P_list_mnta=deepcopy(b[2])
    temp2_freq=deepcopy(b[3])
    temp2_mnta=deepcopy(b[4])
    G_sym=deepcopy(b[5])
    #print 'S_list = ', S_list
    #print 'P_list_freq = ', P_list_freq
    #print 'P_list_mnta = ', P_list_mnta
    #print 'temp2_freq = ', temp2_freq
    #print 'temp2_mnta = ', temp2_mnta
    #print 'G_sym = ', G_sym
    #len_R=0
    #for i in range (0,len(temp2_freq)):
      #if temp2_freq[i]!=None:
        #len_R=len_R+1
      #else:
        #break
 
    for i in range (0,M):  # Pick momenta randomly
      for j in range (0,2):
        k_init[i][j]= (b_h-a_l)*random()+a_l
    #print 'k_init = ', k_init
    poles = Poles(P_list_freq, P_list_mnta,G_sym[0],k_init,t,mu)
    for i in range (0,len(poles)):
      for j in range (0,len(poles[i])):
        if poles[i][j]==None:
          j = len(poles[i])+1
        else:
          for k in range (0,len(poles[i][j])):
            #print poles[i][j][k]
            poles[i][j][k] = f(poles[i][j][k],beta)  # Find the numerical values 

    ########################################
    ################AMI Algebra#############
    ########################################

    res1 = dot_num(S_list[1][0], poles[1][0])
    temp=res1
    for it_m in range (2,M+1):
      o = 0
      for i in range(0, len(S_list[it_m])):
        if S_list[it_m][i] != None:
          o = o + 1
        else:
          i = len(poles[it_m]) + 1
    
      S_new = [None]*(o)
      P_new = [None]*(o)
      for i in range(0, len(S_list[it_m])):
        if S_list[it_m][i] != None:
          S_new[i] = S_list[it_m][i]
          P_new[i] = poles[it_m][i]
  
        else:
          i = len(S_list[it_m]) + 1

      res2 = dot_arr(S_new, P_new)
      temp = cross(res1,res2)
      res1 = temp
  
    G_val = 0.0  # To store numerical value at a given external frequency
    if len(temp2_freq[0])!=0:
      for j in range (0, len(temp)):
        G_val = temp[j]*G_eval(temp2_freq[j], temp2_mnta[j], nu_ext, G_sym[0],k_init,t,mu,beta) + G_val
      #print G_val
    else:
      for j in range (0, len(temp)):
        G_val = temp[j] + G_val  
    out.append(-G_val.real/4**M/pi**(2*M))
    out.append(-G_val.imag/4**M/pi**(2*M))
  return out
      
#End
###################################

###################################
##########Update order#############
###################################
#Begin
def update_order(G,M):
  r=random()
  if r<0.5:
    remove_int_line(G)  # Remove one of the interaction lines
    M=M-1 
    #print G
    #print 'reset ORDER = ', ORDER
    reset_g(G,M)
  if r>0.5:
    add_int_line(G)
    M=M+1
    #print G
    #print 'reset ORDER = ', ORDER
    reset_g(G,M)
  return G,M
#End
###################################

###################################
###########Proposal################
###################################
#Begin
def W(ORDER):
  if ORDER==1:
    return 1.0
  if ORDER==2:
    val=1.0/(b_h-a_l)**(2*ORDER)
    return val
#End
###################################

###################################
#########Acceptence ratio##########
###################################
#Begin
def acp_rt(G_OLD,ORDER_OLD,G_NEW,ORDER_NEW):  # Return the acceptance ratio
  int_old=update_mnta(G_OLD,ORDER_OLD)
  int_old_re=int_old[0]
  int_old_im=int_old[1]
  sign_old_re=numpy.sign(int_old_re)
  sign_old_im=numpy.sign(int_old_im)
  sign_new_re=0
  sign_new_im=0
  if ORDER_NEW>=3 or ORDER_NEW<=0:  # If the order is out of range
    val=0
    return val,sign_old_re,sign_old_im,sign_new_re,sign_new_im
  #print G_OLD
  #print ORDER_OLD
  
  int_new=update_mnta(G_NEW,ORDER_NEW)
  
  int_new_re=int_new[0]
  int_new_im=int_new[1]
  
  sign_new_re=numpy.sign(int_new_re)
  sign_new_im=numpy.sign(int_new_im)
  #print 'int_old_re = ', int_old_re
  #print 'int_new_re = ', int_new_re
  #print 'W = ', W(ORDER_NEW)
  p=(2*min(ORDER_OLD,ORDER_NEW)+1)
  #print 'pB/pA = ', p
  if flag_re==1:
    if int_old_re==0:
      val=1
      return val
    val=abs(int_new_re)/abs(int_old_re)/W(ORDER_NEW)
    if ORDER_NEW>ORDER_OLD:
      val=p*val
    else:
      val=val/p
  return val,sign_old_re,sign_old_im,sign_new_re,sign_new_im
#End
###################################
     
###################################
##########Update function##########
###################################
#Begin
def new_state(G_old,ORDER_old):  # Accept or reject the proposed state
  OLD_diag=G_old.copy()
  proposed_diag_ord=update_order(G_old,ORDER_old)
  diag_proposed=proposed_diag_ord[0]
  ORDER_proposed=proposed_diag_ord[1]
  #print 'g_old = ',  G_old
  #print 'm_old = ',  ORDER_old
  #print 'g_proposed = ',  diag_proposed
  #print 'm_proposed = ',  ORDER_proposed
  R_signs=acp_rt(OLD_diag,ORDER_old,diag_proposed,ORDER_proposed)
  #print R_signs 
  R=R_signs[0]
  #print 'R = ', R
  r=random()
  #print 'r = ', r
  if R <= r:
    NEW_G = OLD_diag
    NEW_ORDER=ORDER_old
    NEW_sign_re=R_signs[1]
    NEW_sign_im=R_signs[2]
  else:  
    NEW_G=diag_proposed
    NEW_ORDER=ORDER_proposed 
    NEW_sign_re=R_signs[3]
    NEW_sign_im=R_signs[4]
  if flag_re==1:
    sign=NEW_sign_re
  if flag_re==0:
    sign=NEW_sign_im  
  #print NEW_G,NEW_ORDER,sign 
  return NEW_G,NEW_ORDER,sign
  
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
next_state=[]
flag_re = 1  # Global variable
a_l=-pi  # Lower limit of the integrals
b_h=pi   # Upper limit of the integrals
m_old=2   # Order of the initial state
old_diag=generate_g_2(m_old)  # generate initial state
#old_diag=g1.copy()
print 'old_diag = ', old_diag

for i in range (0,1000):
  next_state.append(new_state(old_diag,m_old))
  #print 'next_diag = ', next_state[0][0]
  #print 'next_order = ', next_state[0][1]
  #print next_state[0][2]
  reset_g(next_state[0][0],next_state[0][1])
  next_diag=next_state[0][0].copy()
  next_order=next_state[0][1]
  next_sign=next_state[0][2]
  old_diag=next_state[0][0].copy()
  m_old=next_state[0][1]
  #print 'old_diag = ', old_diag
  print 'm_old = ', m_old
  reset_g(old_diag,m_old)  
  next_state=[]
  
  


























#CHECK = checks(new_g[0],new_g[1])
#flag_AMI = CHECK[0]
#print new_g[1]
#print update_mnta(new_g[0],new_g[1])
'''
if flag_AMI==1:
  
  b= AMI_arrays_out(CHECK[1],new_g[1])   # Generate S, P, and R arrays and store in b
  print("--- Construction_time ---", (time.time() - start_time))  # Compute the construction time
  print 'S = ', b[0]
  print 'P_freq = ', b[1]
  print 'P_mnta = ', b[2]
  print 'R_freq = ', b[3]
  print 'R_mnta = ', b[4]
'''



