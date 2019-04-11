##################################
##########Import Libraries########
##################################
#Begin
from igraph import *
from random import *
import numpy
#End
###################################

###################################
#######External Frequencies########
###################################
#Begin
def EXTERNAL_F(M):
  ext_f = [0]*(M+1+1) 
  ext_f[M+1] = 1
  return ext_f
#End
####################################

###################################
#####Independent Frequencies#######
###################################
#Begin
def INTERNAL_F(M):
  ind_f=[[0]*(M+1+1)for i in range (0,M+1)]
  for j in range (0,M+1):
    ind_f[j][j]=1
  return ind_f
#End
####################################

###################################
###Counting Non-labelled edges#####
###################################
#Begin
def cnt_n_l(G,M,V):
  cnt=0               # To store number of the nonlabelled edges for a specific node
  #print "node = ", node           # Check node
  for i in range (0,len(G.adjacent(V, mode=IN))):  # Compute number of nonlabelled ingoing edges	
    if G.es["Label"][G.adjacent(V, mode=IN)[i]]==[None]*(M+1+1):	
      cnt=cnt+1
  for i in range (0,len(G.adjacent(V, mode=OUT))):  # Compute number of nonlabelled outgoing edges	
    if G.es["Label"][G.adjacent(V, mode=OUT)[i]]==[None]*(M+1+1):	
      cnt=cnt+1
  return cnt
#End
####################################

###################################
####Labelling for non-count=2######
###################################
#Begin
def label_2(G,M,V,CNT,NS): # G: graph, M: order, V: the node, CNT: number of the nonlabelled edges connected to V, NS: independent frequency counter
  ind_f=INTERNAL_F(M)
  #print 'ind_f =', ind_f
  Ns= NS
  Cnt = CNT
  TF=False
  for i in range(0,len(G.adjacent(V, mode=IN))):  # For ingoing edges
    if Cnt>1:
      if G.es[G.adjacent(V, mode=IN)[i]]["Label"]==[None]*(M+1+1):
        if (G.es[G.adjacent(V, mode=IN)[i]]["F_or_B"]!=0) and (Ns<=M+1-1): # Assign an independent frequency to a fermionic line
          #label[G.adjacent(V, mode=IN)[i]]=ind_f[Ns]
          G.es[G.adjacent(V, mode=IN)[i]]["Label"]=ind_f[Ns]
          #nus[Ns]=ind_f[Ns]
          Ns = Ns+1
          Cnt = Cnt-1  
          TF=True
  for i in range(0,len(G.adjacent(V, mode=OUT))):  # For ingoing edges
    if Cnt>1:
      if G.es[G.adjacent(V, mode=OUT)[i]]["Label"]==[None]*(M+1+1):
        if (G.es[G.adjacent(V, mode=OUT)[i]]["F_or_B"]!=0) and (Ns<=M+1-1): # Assign an independent frequency to a fermionic line
          #label[G.adjacent(V, mode=OUT)[i]]=ind_f[Ns]
          G.es[G.adjacent(V, mode=OUT)[i]]["Label"]=ind_f[Ns]
          #nus[Ns]=ind_f[Ns]
          Ns = Ns+1
          Cnt = Cnt-1  
          TF=True
  return G, Cnt, Ns   
#End
###################################

###################################
####Labelling for non-count=1######
###################################
#Begin
def label_1(G,M,V):
      index=None
      VAL = [0]*(M+1+1)  # Auxiliary variable to store a label temporarily
      for i in range(0,len(G.adjacent(V, mode=IN))):  
        if (G.es[G.adjacent(V, mode=IN)[i]]["Label"]==[None]*(M+1+1)):  # Find the unlabelled ingoing edge 
          index=G.adjacent(V, mode=IN)[i]
          G.es[G.adjacent(V, mode=IN)[i]]["Label"]=[0]*(M+1+1)
          sign=-1  # To consider flow direction
    
        for k in range(0,M+1+1): 
          VAL[k] = VAL[k] + G.es[G.adjacent(V, mode=IN)[i]]["Label"][k] # Compute total ingoing flow

      for i in range(0,len(G.adjacent(V, mode=OUT))): 
        if (G.es[G.adjacent(V, mode=OUT)[i]]["Label"]==[None]*(M+1+1)):  # Find the unlabelled outgoing edge
          G.es[G.adjacent(V, mode=OUT)[i]]["Label"]=[0]*(M+1+1)
          index=G.adjacent(V, mode=OUT)[i]
          sign=1  # To consider flow direction
        for k in range(0,M+1+1):
          VAL[k] = VAL[k] - G.es[G.adjacent(V, mode=OUT)[i]]["Label"][k] # Compute total outgoing flow
      if index!=None:
        for k in range(0,M+1+1): # Compute the last label
          VAL[k]=sign*VAL[k]
        G.es[index]["Label"] = VAL # Assign the label
      
      return G      
#End
####################################   


####################################
##########Array Comparison##########
####################################
#Begin
def comp(a,b):
  oh = True
  for i in range (0,len(a)):
    if a[i]!=b[i]:
      oh = False
      break
  return oh
#End
#####################################

####################################
###########Reset Labels#############
####################################
def reset_g(G,M):
  #Begin
  l = L(M)
  EXT_F=EXTERNAL_F(M)
  G.es[0]["Label"]=EXT_F
  for i in range(1,l):
    if G.es[i]["INT_or_EXT"]==0:
      G.es[i]["Label"]=EXT_F
  return G
  #End
####################################

###################################
####Labelling Procedure: Random####
###################################
#Begin
def label_ran(G,M,MAX_V):
  n_v = N_V(M)
  l = L(M)
  G=reset_g(G,M)   # Reset g
  lbl_num=0
  ns=0
  try_label=1 
  count_non_label=0

  while lbl_num<l:   # Until all the lines are labelled
  
    for i in range(0,MAX_V):  # Pick max_v nodes
      node = randint(0, n_v-1)            # Picking a random node
      count_non_label = cnt_n_l(G,M,node)  # Computing the nonlabelled edges of the node
      #print 'node ', node, 'count ', count_non_label  
      #print 'count_non_label = ', count_non_label
  
      if (count_non_label==2):   # Do the calculations for nodes with 2 unnkown edges
        temp=label_2(G,M,node,count_non_label,ns) 
        G=temp[0]
        count_non_label=temp[1]
        #if temp[2]:
        ns = temp[2]

      if count_non_label==1:  # If there is only one unknown edge find the last one considering conservation laws
        G=label_1(G,M,node)
      check = 0
      while check == 0:  # Label all the nodes with only one unknown edge
        check=1
        for node in range (0,n_v):
          count_non_label = cnt_n_l(G,M,node) 
          if (count_non_label==1):
            G=label_1(G,M,node)  
            check=0
      
      ######################################
      #####Number of the labelled edges#####
      ######################################
      lbl_num=0
      for i in range(0,l):
        if (comp(G.es[i]["Label"],[None]*l)==False):
          lbl_num=lbl_num+1    
          #print 'lbl_num =', lbl_num
          #print G.es["Label"]
          #End
      if lbl_num==l:
        break
    #######################################

    ####################################
    ############Reset Labels############
    ####################################
    #Begin
    if lbl_num<l:
      G.es["Label"] = [[None]*(M+1+1)]*l  # Define "None" labels 
      G=reset_g(G,M)
      lbl_num=0
      ns=0
      try_label = try_label+1
      #print g.es["Label"]
    #End
    ####################################

  #End
  return G, try_label
#End
#####################################

###################################
##Labelling Procedure: Systematic##
###################################
#Begin
def label_sys(G,M):
  lbl_num=0    # Numebr of labelled edges
  ns=0         # Number of the assigned independent frequencies
  counter=0    # Number of the picked nodes in the path 
  out = [None]*(M+1)   # Nodes list
  n_v = N_V(M)   # Total number of the nodes
  l = L(M)   # Total number of the lines 
  count_non_label=0  # Number of the nonlabelled edges adjacent to the node
  G=reset_g(G,M)   # Reset g
  for i in range(0,n_v):
    count_non_label = cnt_n_l(G,M,i)  # Computing the nonlabelled edges of the node
    if count_non_label==2:
      node = i
      #if randint(0,1)==0:
      break

  while lbl_num<l:   # Until all the lines are labelled
    if counter==M+1:  # If labelling was not successful do it again
      G=generate(M)   # Generate g
      G=reset_g(G,M)  # Reset labels, and assign external lines' labels
      lbl_num=0
      ns=0
      counter=0 
      for i in range(0,n_v):
        count_non_label = cnt_n_l(G,M,i)  # Computing the nonlabelled edges of the node
        if count_non_label==2:
          node = i
          break 
    G.vs[node]["visited"]=1  # The node is visited
    count_non_label = cnt_n_l(G,M,node)  # Computing the nonlabelled edges of the node
    
    out[counter] = node  # Put the node in the list
    counter = counter + 1  
    print node
  
    if (count_non_label==2):   # Do the calculations for nodes with 2 unnkown edges
      temp=label_2(G,M,node,count_non_label,ns) 
      G=temp[0]
      count_non_label=temp[1]
      ns = temp[2]

    if count_non_label==1:  # If there is only one unknown edge find the last 
                            # one considering conservation laws
      G=label_1(G,M,node)
    check = 0
    while check == 0:  # Label all the nodes with only one unknown edge
      check=1
      for j in range (0,n_v):  # Sweep all the nodes
        count_non_label = cnt_n_l(G,M,j) 
        if (count_non_label==1):
          G=label_1(G,M,j)  
          check=0
      
    ######################################
    #####Number of the labelled edges#####
    ######################################
    lbl_num=0
    for i in range(0,l):  # Sweep all the edges
      if (comp(G.es[i]["Label"],[None]*l)==False):
        lbl_num=lbl_num+1    

    if lbl_num==l:   # If the labelling is completed
      break
    #######################################
    
    ######################################
    ##########Pick the next node##########
    ######################################
    neighbor_non_visit=[None]*2
    cnt_non_visit=0
    for i in range(0,len(G.neighbors(node))):
      if (G.vs[G.neighbors(node)[i]]["visited"]==0) and (cnt_n_l(G,M,G.neighbors(node)[i])==2):
        neighbor_non_visit[cnt_non_visit]=G.neighbors(node)[i]
        cnt_non_visit = cnt_non_visit+1       
    if cnt_non_visit==1:
      node=neighbor_non_visit[0]
    if cnt_non_visit==2:
      node=neighbor_non_visit[randint(0,1)]
    if cnt_non_visit==0:
      for i in range(0,n_v):  # Pick next node by sweeping the nodes
      #while True:  # Pick nodes randomly
        #i = randint(0,n_v-1)
        count_non_label = cnt_n_l(G,M,i)  # Computing the nonlabelled edges of the node
        if count_non_label==2 and G.vs[i]["visited"]==0:
          node = i
          break   
  #End  
  ################################# 
  return G, out
#End
#####################################

###################################
########Labelling Procedure########
###################################
#Begin
def label_f_line(G,M):
  lbl_num=0
  ns=0
  counter=0 
  out = [None]*(M+1)
  n_v = N_V(M)
  l = L(M)
  count_non_label=0
  G=reset_g(G,M)   # Reset g
  for i in range(0,n_v):  # Pick the first node
    count_non_label = cnt_n_l(G,M,i)  # Computing the nonlabelled edges of the node
    if count_non_label==2:
      node = i
      break
  
  while lbl_num<l:   # Until all the lines are labelled
    for i in range(0,n_v):  # Sweep loop
      #print 'node = ', node
      #print 'counter = ', counter
      G.vs[node]["visited"]=1   # The node is visited
      count_non_label = cnt_n_l(G,M,node)  # Computing the nonlabelled edges of the node
      #print node 
      out[counter] = node  # One of the nodes in the path
      #print 'count_non_label = ', count_non_label
  
      if (count_non_label==2):   # Do the calculations for nodes with 2 unnkown edges
        temp=label_2(G,M,node,count_non_label,ns) 
        G=temp[0]
        count_non_label=temp[1]
        ns = temp[2]
      
        counter = counter + 1  # Go to the next node index

      if count_non_label==1:  # If there is only one unknown edge find the last one considering conservation laws
        G=label_1(G,M,node)
      check = 0
      while check == 0:  # Label all the nodes with only one unknown edge
        check=1
        for j in range (0,n_v):
          count_non_label = cnt_n_l(G,M,j) 
          if (count_non_label==1):
            G=label_1(G,M,j)  
            check=0
      
      ######################################
      #####Number of the labelled edges#####
      ######################################
      lbl_num=0
      for i in range(0,l):
        if (comp(G.es[i]["Label"],[None]*l)==False):
          lbl_num=lbl_num+1    
      if lbl_num==l:  # If the labelling is complete
        break

      ######################################
      #######Pick the next node in loop#####
      ######################################
      x = G.neighbors(node, mode="OUT")  # Out neighbors
      #print x
      flag = True  # Auxiliary variable to choose the fermionic line
      for i in range (0,len(G.neighbors(node, mode="OUT"))):
        if G.es[G.get_eid(node,x[i])]["F_or_B"] == 1:  # If the outgoing line is fermionic
          node = G.neighbors(node, mode="OUT")[i] # Go to the next node
          flag = False 
          break
      if flag: # If the next has not been determined
        x = G.neighbors(node, mode="IN")
        for i in range (0,len(G.neighbors(node, mode="IN"))):
          if G.es[G.get_eid(x[i],node)]["F_or_B"] == 1:  # If the ingoing line is fermionic
            node = G.neighbors(node, mode="IN")[i]  # Go to the next node
            break      
    
        
    #End  
    ##############################
    
    ########################################
    ####Go to the next loop if necessary####
    ########################################
    
    if lbl_num<l:  # If the diagram has more than one loop
      #print "Go to the next loop" 
      for i in range(0,n_v):
        if G.vs[i]["visited"]==0:
          node = i
  return G, out
#End
#####################################


##############################
####Compute l_f, l_i,l,n_v####
##############################

def L_F(M):
  return 2*(M+1)

def L_I(M):
  return M

def L(M):
  return L_F(M)+L_I(M)+2

def N_V(M):
  return 2*(M+2)

#End
################################

#############################
##Define a 0th order Graph###
#############################
#Begin
def generate_g_0(M):  # M: order of diagram
  n_v = N_V(M)
  l = L(M)
  G = Graph(directed=True)
  G.add_vertices(n_v) 

  G.add_edges([(0,1)])
  G.add_edges([(1,2)])
  G.add_edges([(2,1)])
  G.add_edges([(2,3)])

  #####################################
  ############Attributes###############
  #####################################
  #Begin
  G.vs["name1"] = ["0", "1", "2", "3"]
  G.vs["visited"] = [1,0,0,1]
  G.es["name2"] = ["0-B", "1-F", "2-F", "3-B"]
  G.es["F_or_B"] = [0,1,1,0]
  G.es["INT_or_EXT"] = [0,1,1,0]
  G.es["Label"] = [[None]*(M+1+1)]*l  # Define "None" labels 
  #End
  #####################################
  
  return G

#End
#######################################

#############################
##Define a 1st order Graph###
#############################
#Begin
def generate_g_1(M):  # M: order of diagram
  n_v = N_V(M)
  l = L(M)
  G = Graph(directed=True)
  G.add_vertices(n_v) 

  G.add_edges([(0,1)])
  G.add_edges([(1,2)])
  G.add_edges([(4,1)])
  G.add_edges([(2,3)])
  G.add_edges([(3,2)])
  G.add_edges([(3,4)])
  G.add_edges([(4,5)])

  #####################################
  ############Attributes###############
  #####################################
  #Begin
  G.vs["name1"] = ["0", "1", "2", "3", "4", "5"]
  G.vs["visited"] = [1,0,0,0,0,1]
  G.es["name2"] = ["0-B", "1-F", "2-F", "3-F", "4-B", "5-F", "6-B"]
  G.es["F_or_B"] = [0,1,1,1,0,1,0]
  G.es["INT_or_EXT"] = [0,1,1,1,1,1,0]
  G.es["Label"] = [[None]*(M+1+1)]*l  # Define "None" labels 
  #End
  #####################################
  
  return G

#End
######################################


#############################
##Define a 2nd order Graph###
#############################
#Begin
def generate_g_2_old(M):  # M: order of diagram
  n_v = N_V(M)
  l = L(M)
  G = Graph(directed=True)
  G.add_vertices(n_v) 

  G.add_edges([(0,1)])
  G.add_edges([(1,2)])
  G.add_edges([(6,1)])
  G.add_edges([(2,3)])
  G.add_edges([(2,4)])
  G.add_edges([(3,4)])
  G.add_edges([(3,6)])
  G.add_edges([(4,5)])
  G.add_edges([(5,6)])
  G.add_edges([(5,7)])

  #####################################
  ############Attributes###############
  #####################################
  #Begin
  G.vs["name1"] = ["0", "1", "2", "3", "4", "5", "6", "7"]
  G.vs["visited"] = [1,0,0,0,0,0,0,1]
  G.es["name2"] = ["0-B", "1-F", "2-F", "3-F", "4-B", "5-F", "6-B", "7-F", "8-F", "9-B" ]
  G.es["INT_or_EXT"] = [0,1,1,1,1,1,1,1,1,0]
  G.es["F_or_B"] = [0,1,1,1,0,1,0,1,1,0]
  G.es["Label"] = [[None]*(M+1+1)]*l  # Define "None" labels 
  #End
  #####################################
  
  return G

#End
######################################

#############################
##Define a 2nd order Graph###
#############################
#Begin
def generate_g_2(M):  # M: order of diagram
  #Begin
  G = Graph(directed=True)
  n_v = N_V(M)
  l = L(M)
  G.add_vertices(n_v) 

  G.add_edges([(0,2)])  #B
  G.add_edges([(2,3)])  #F
  G.add_edges([(3,4)])  #B
  G.add_edges([(4,5)])  #F
  G.add_edges([(5,4)])  #F
  G.add_edges([(5,6)])  #B
  G.add_edges([(6,2)])  #F
  G.add_edges([(3,7)])  #F
  G.add_edges([(7,6)])  #F
  G.add_edges([(7,1)])  #B
  
  #End
  #####################################

  #####################################
  ############Attributes###############
  #####################################
  #Begin
  G.vs["name1"] = ["0", "1", "2", "3", "4", "5", "6", "7"]
  G.vs["visited"] = [1,1,0,0,0,0,0,0]
  G.es["name2"] = ["0-B", "1-F", "2-B", "3-F", "4-F", "5-B", "6-F", "7-F", "8-F", "9-B"]
  # Bosonic line: 0, # Femionic line: 1
  G.es["F_or_B"] = [0,1,0,1,1,0,1,1,1,0]
  # External line: 0, # Internal line: 1
  G.es["INT_or_EXT"] = [0,1,1,1,1,1,1,1,1,0]
  G.es["Label"] = [[None]*(M+1+1)]*l  # Define "None" labels 
  #End
  #####################################
  
  return G
#####################################

#############################
##Define a 3rd order Graph###
#############################
#Begin
def generate_g_3(M):  # M: order of diagram
  #Begin
  G = Graph(directed=True)
  n_v = N_V(M)
  l = L(M)
  G.add_vertices(n_v) 

  G.add_edges([(0,1)])
  G.add_edges([(1,2)])
  G.add_edges([(8,1)])
  G.add_edges([(2,3)])
  G.add_edges([(2,4)])
  G.add_edges([(3,4)])
  G.add_edges([(3,6)])
  G.add_edges([(4,5)])
  G.add_edges([(5,6)])
  G.add_edges([(5,7)])
  G.add_edges([(6,7)])
  G.add_edges([(7,8)])
  G.add_edges([(8,9)])

  #End
  #####################################

  #####################################
  ############Attributes###############
  #####################################
  #Begin
  G.vs["name1"] = ["0", "1", "2", "3", "4", "5", "6", "7","8","9"]
  G.vs["visited"] = [1,0,0,0,0,0,0,0,0,1]
  G.es["name2"] = ["0-B", "1-F", "2-F", "3-F", "4-B", "5-F", "6-B", "7-F", "8-F", "9-B","10-F", "11-F", "12-B"]
  # Bosonic line: 0, # Femionic line: 1
  G.es["F_or_B"] = [0,1,1,1,0,1,0,1,1,0,1,1,0]
  # External line: 0, # Internal line: 1
  G.es["INT_or_EXT"] = [0,1,1,1,1,1,1,1,1,1,1,1,0]
  G.es["Label"] = [[None]*(M+1+1)]*l  # Define "None" labels 
  #End
  #####################################
  
  return G
#####################################

#############################
##Define a 4th order Graph###
#############################
#Begin
def generate_g_4(M):  # M: order of diagram
  #Begin
  G = Graph(directed=True)
  n_v = N_V(M)
  l = L(M)
  G.add_vertices(n_v) 

  G.add_edges([(0,2)])
  G.add_edges([(2,3)])
  G.add_edges([(11,2)])
  G.add_edges([(3,4)])
  G.add_edges([(3,8)])
  G.add_edges([(4,5)])
  G.add_edges([(4,10)])
  G.add_edges([(5,6)])
  G.add_edges([(6,5)])
  G.add_edges([(6,7)])
  G.add_edges([(7,8)])
  G.add_edges([(7,11)])
  G.add_edges([(8,9)])
  G.add_edges([(9,1)])
  G.add_edges([(9,10)])
  G.add_edges([(10,11)])

  #End
  #####################################

  #####################################
  ############Attributes###############
  #####################################
  #Begin
  G.vs["name1"] = ["0", "1", "2", "3", "4", "5", "6", "7","8","9","10","11"]
  G.vs["visited"] = [1,1,0,0,0,0,0,0,0,0,0,0]
  G.es["name2"] = ["0-B", "1-F", "2-F", "3-F", "4-B", "5-F", "6-B", "7-F", "8-B", "9-F","10-F", "11-B", "12-F","13-B","14-F","15-F"]
  # Bosonic line: 0, # Femionic line: 1
  G.es["F_or_B"] = [0,1,1,1,0,1,0,1,0,1,1,0,1,0,1,1]
  # External line: 0, # Internal line: 1
  G.es["INT_or_EXT"] = [0,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1]
  G.es["Label"] = [[None]*(M+1+1)]*l  # Define "None" labels 
  #End
  #####################################
  
  return G
#####################################


#############################
##Define a 5th order Graph###
#############################
#Begin
def generate_g_5(M):  # M: order of diagram
  #Begin
  G = Graph(directed=True)
  n_v = N_V(M)
  l = L(M)
  G.add_vertices(n_v) 

  G.add_edges([(0,2)])
  G.add_edges([(2,3)])
  G.add_edges([(13,2)])
  G.add_edges([(3,4)])
  G.add_edges([(3,6)])
  G.add_edges([(4,5)])
  G.add_edges([(4,7)])
  G.add_edges([(5,6)])
  G.add_edges([(5,11)])
  G.add_edges([(6,7)])
  G.add_edges([(7,8)])
  G.add_edges([(8,9)])
  G.add_edges([(8,10)]) 
  G.add_edges([(9,1)])
  G.add_edges([(9,10)])
  G.add_edges([(10,11)])
  G.add_edges([(11,12)])
  G.add_edges([(12,13)])
  G.add_edges([(13,12)])
  
  
  
  #End
  #####################################

  #####################################
  ############Attributes###############
  #####################################
  #Begin
  G.vs["name1"] = ["0", "1", "2", "3", "4", "5", "6", "7","8","9","10","11","12","13"]
  G.vs["visited"] = [1,1,0,0,0,0,0,0,0,0,0,0,0,0]
  G.es["name2"] = ["0-B", "1-F", "2-F", "3-F", "4-B", "5-F", "6-B", "7-F", "8-B", "9-F","10-F", "11-F", "12-B","13-B","14-F","15-F","16-F","17-F","18-B"]
  # Bosonic line: 0, # Femionic line: 1
  G.es["F_or_B"] = [0,1,1,1,0,1,0,1,0,1,1,1,0,0,1,1,1,1,0]
  # External line: 0, # Internal line: 1
  G.es["INT_or_EXT"] = [0,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1]
  G.es["Label"] = [[None]*(M+1+1)]*l  # Define "None" labels 
  #End
  #####################################
  
  return G
#####################################

#############################
##Define a 6th order Graph###
#############################
#Begin
def generate_g_6(M):  # M: order of diagram
  #Begin
  G = Graph(directed=True)
  n_v = N_V(M)
  l = L(M)
  G.add_vertices(n_v) 

  G.add_edges([(0,2)])
  G.add_edges([(2,3)])
  G.add_edges([(15,2)])
  G.add_edges([(3,4)])
  G.add_edges([(3,5)])
  G.add_edges([(4,5)])
  G.add_edges([(4,6)])
  G.add_edges([(5,6)])
  G.add_edges([(6,7)])
  G.add_edges([(7,8)])
  G.add_edges([(7,9)])
  G.add_edges([(8,9)])
  G.add_edges([(8,12)]) 
  G.add_edges([(9,10)])
  G.add_edges([(10,1)])
  G.add_edges([(10,11)])
  G.add_edges([(11,12)])
  G.add_edges([(11,13)])
  G.add_edges([(12,13)])
  G.add_edges([(13,14)])
  G.add_edges([(14,15)])
  G.add_edges([(15,14)])
  
  #End
  #####################################

  #####################################
  ############Attributes###############
  #####################################
  #Begin
  G.vs["name1"] = ["0", "1", "2", "3", "4", "5", "6", "7","8","9","10","11","12","13","14","15"]
  G.vs["visited"] = [1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
  G.es["name2"] = ["0-B", "1-F", "2-F", "3-F", "4-B", "5-F", "6-B", "7-F", "8-F", "9-F","10-B", "11-F", "12-B","13-F","14-B","15-F","16-F","17-B","18-F","19-F","20-F","21-B"]
  # Bosonic line: 0, # Femionic line: 1
  G.es["F_or_B"] = [0,1,1,1,0,1,0,1,1,1,0,1,0,1,0,1,1,0,1,1,1,0]
  # External line: 0, # Internal line: 1
  G.es["INT_or_EXT"] = [0,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1]
  G.es["Label"] = [[None]*(M+1+1)]*l  # Define "None" labels 
  #End
  #####################################
  
  return G
#####################################

#############################
##Define a 7th order Graph###
#############################
#Begin
def generate_g_7(M):  # M: order of diagram
  G = Graph(directed=True)
  n_v = N_V(M)
  l = L(M)
  G.add_vertices(n_v) 
  G.add_edges([(0,1)])
  G.add_edges([(1,2)])
  G.add_edges([(15,1)])
  G.add_edges([(2,3)])
  G.add_edges([(2,4)])
  G.add_edges([(3,4)])
  G.add_edges([(3,12)])
  G.add_edges([(4,5)])
  G.add_edges([(5,6)])
  G.add_edges([(5,7)])
  G.add_edges([(6,7)])
  G.add_edges([(6,9)])
  G.add_edges([(7,8)])
  G.add_edges([(8,9)])
  G.add_edges([(8,14)])
  G.add_edges([(9,10)])
  G.add_edges([(10,11)])
  G.add_edges([(10,16)])
  G.add_edges([(16,11)])
  G.add_edges([(16,17)])
  G.add_edges([(11,12)])
  G.add_edges([(12,13)])
  G.add_edges([(13,14)])
  G.add_edges([(13,15)])
  G.add_edges([(14,15)])
  
  #End
  #####################################
  
  #####################################
  ############Attributes###############
  #####################################
  #Begin
  G.vs["name1"] = ["0", "1", "2", "3", "4", "5", "6", "7","8","9","10","11","12","13","14","15","16","17"]
  G.vs["visited"] = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1]
  G.es["name2"] = ["0-B", "1-F", "2-F", "3-F", "4-B", "5-F", "6-B", "7-F", "8-F", "9-B","10-F", "11-B", "12-F", "13-F", "14-B", "15-F", "16-B", "17-F","18-F","19-B","20-F","21-F","22-F","23-B","24-F"]
  # Bosonic line: 0, # Femionic line: 1
  G.es["F_or_B"] = [0,1,1,1,0,1,0,1,1,0,1,0,1,1,0,1,0,1,1,0,1,1,1,0,1]
  # External line: 0, # Internal line: 1
  G.es["INT_or_EXT"] = [0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1]
  G.es["Label"] = [[None]*(M+1+1)]*l  # Define "None" labels 
  #End
  #####################################
  
  return G

#End
########################################

#############################
##Define a Graph of order M##
#############################
def generate(M):
#Begin
  if M==0:
    return generate_g_0(M)

  if M==1:
    return generate_g_1(M)

  if M==2:
    return generate_g_2(M)

  if M==3:
    return generate_g_3(M)

  if M==4:
    return generate_g_4(M)

  if M==5:
    return generate_g_5(M)

  if M==6:
    return generate_g_6(M)

  if M==7:
    return generate_g_7(M)

#End
#############################

