##################################
##########Import Libraries########
##################################
#Begin
from random import *
from igraph import *
import numpy
#End
###################################

##############################
####Compute l_f, l_i,l,n_v####
##############################

def L_F(M):
  return 2*M+1

def L_I(M):
  return M

def L(M):
  return L_F(M)+L_I(M)

def N_V(M):
  return 2*M+2

#End
################################

########################################
##Define a 2nd order self-energy graph##
########################################
#Begin
def generate_g_1_con(M):  # M: order of diagram
  n_v = N_V(M)
  l = L(M)
  G = Graph(directed=True)
  G.add_vertices(n_v) 

  G.add_edges([(2,3)])
  G.add_edges([(3,1)])
  G.add_edges([(0,3)])
  G.add_edges([(2,2)])

  #####################################
  ############Attributes###############
  #####################################
  #Begin
  G.vs["name1"] = [0, 1, 2, 3]
  G.vs["visited"] = [1,1,0,0]
  G.vs["Spin"] = [None]*n_v
  G.es["name2"] = ["0-F", "1-F", "2-B", "3-F"]
  G.es["F_or_B"] = [0,1,1,1]
  G.es["INT_or_EXT"] = [1,0,0,1]
  G.es["Label"] = [[None]*(M+1)]*l  # Define "None" labels
   
  #End
  #####################################
  
  return G

#End
#############################################

########################################
##Define a 2nd order self-energy graph##
########################################
#Begin
def generate_g_2(M):  # M: order of diagram
  n_v = N_V(M)
  l = L(M)
  G = Graph(directed=True)
  G.add_vertices(n_v) 

  G.add_edges([(0,2)])
  G.add_edges([(2,3)])
  G.add_edges([(3,4)])
  G.add_edges([(4,3)])
  G.add_edges([(4,5)])
  G.add_edges([(2,5)])
  G.add_edges([(5,1)])

  #####################################
  ############Attributes###############
  #####################################
  #Begin
  G.vs["name1"] = [0, 1, 2, 3, 4, 5]
  G.vs["visited"] = [1,1,0,0,0,0]
  G.vs["Spin"] = [None]*n_v
  #G.es["name2"] = ["0-F", "1-B", "2-F", "3-F", "4-B", "5-F", "6-F"]
  G.es["F_or_B"] = [1,0,1,1,0,1,1]
  G.es["INT_or_EXT"] = [0,1,1,1,1,1,0]
  G.es["Label"] = [[None]*(M+1)]*l  # Define "None" labels
   
  #End
  #####################################
  
  return G

#End
######################################

###################################
###########Add a tadpole###########
###################################
#Begin
def add_tad(G,M):
  f_list=[]
  NV=2*M+2  # number of vertices
  l=3*M+1  # number of lines
  for i in range (0,l):   # Sweep all the lines
    if G.es[i]['F_or_B'] == 1 and G.es[i]["INT_or_EXT"]==1:   # If the line is fermionic and internal
      f_list.append(i) 
  #print 'f_list = ', f_list
  lif=len(f_list)
  f = randint(0,lif-1)  # Choose on of the internal solid lines
  #print f
  V1 = G.get_edgelist()[f_list[f]][0]  # The first node connected to f1
  print 'V1 = ', V1
  V2 = G.get_edgelist()[f_list[f]][1]  # The second node connected to f1
  print 'V2 = ', V2
  G.add_vertices(1)
  G.vs[NV]['name1'] = NV  # Update 'name1' attribute of the new node
  G.vs[NV]["visited"] = 0       # Update 'visited' attribute of the new node
  #print 'Diagram after one node = ', G
  #print 'name1 after one node = ', G.vs['name1']
  #print 'visited after one node = ', G.vs['visited']
  #print 'F_or_B after one node = ', G.es['F_or_B']
  #print 'INT_or_EXT after one node = ', G.es["INT_or_EXT"]
  G.delete_edges(f_list[f])  # Delete f
  #print 'Diagram after deleting one edge = ', G
  #print 'name1 after deleting one edge = ', G.vs['name1']
  #print 'visited after deleting one edge = ', G.vs['visited']
  #print 'F_or_B after deleting one edge = ', G.es['F_or_B']
  #print 'INT_or_EXT after deleting one edge = ', G.es["INT_or_EXT"]
  G.add_edges([(V1,NV)])     # Add first edge between V1 and new node
  G.es[G.get_eid(V1,NV)]['F_or_B'] = 1   # Update 'F_or_B' attribute of the new edge
  G.es[G.get_eid(V1,NV)]['INT_or_EXT'] = 1   # Update 'INT_or_EXT' attribute of the new edge
  #print 'Diagram after adding first edge = ', G
  #print 'name1 after adding first edge = ', G.vs['name1']
  #print 'visited after adding first edge = ', G.vs['visited']
  #print 'F_or_B after adding first edge = ', G.es['F_or_B']
  #print 'INT_or_EXT after adding first edge = ', G.es["INT_or_EXT"]
  G.add_edges([(NV,V2)])     # Add second edge between new node and V2
  G.es[G.get_eid(NV,V2)]['F_or_B'] = 1   # Update 'F_or_B' attribute of the new edge
  G.es[G.get_eid(NV,V2)]['INT_or_EXT'] = 1   # Update 'INT_or_EXT' attribute of the new edge
  #print 'Diagram after adding second edge = ', G
  #print 'name1 after adding second edge = ', G.vs['name1']
  #print 'visited after adding second edge = ', G.vs['visited']
  #print 'F_or_B after adding second edge = ', G.es['F_or_B']
  #print 'INT_or_EXT after adding second edge = ', G.es["INT_or_EXT"]
  G.add_vertices(1)
  G.vs[NV+1]['name1'] = NV+1  # Update 'name1' attribute of the new node
  G.vs[NV+1]["visited"] = 0       # Update 'visited' attribute of the new node
  #print 'Diagram after second node = ', G
  #print 'name1 after one node = ', G.vs['name1']
  #print 'visited after one node = ', G.vs['visited']
  #print 'F_or_B after one node = ', G.es['F_or_B']
  #print 'INT_or_EXT after one node = ', G.es["INT_or_EXT"]
  G.add_edges([(NV,NV+1)])     # Add first edge between V1 and new node
  G.es[G.get_eid(NV,NV+1)]['F_or_B'] = 0   # Update 'F_or_B' attribute of the new edge
  G.es[G.get_eid(NV,NV+1)]['INT_or_EXT'] = 1   # Update 'INT_or_EXT' attribute of the new
  #print 'Diagram after adding third edge = ', G
  #print 'name1 after adding third edge = ', G.vs['name1']
  #print 'visited after adding third edge = ', G.vs['visited']
  #print 'F_or_B after adding third edge = ', G.es['F_or_B']
  #print 'INT_or_EXT after adding third edge = ', G.es["INT_or_EXT"]
  G.add_edges([(NV+1,NV+1)])     # Add first edge between V1 and new node
  G.es[G.get_eid(NV+1,NV+1)]['F_or_B'] = 1   # Update 'F_or_B' attribute of the new edge
  G.es[G.get_eid(NV+1,NV+1)]['INT_or_EXT'] = 1   # Update 'INT_or_EXT' attribute of the new
  #print 'Diagram after adding fourth edge = ', G
  #print 'name1 after adding fourth edge = ', G.vs['name1']
  #print 'visited after adding fourth edge = ', G.vs['visited']
  #print 'F_or_B after adding fourth edge = ', G.es['F_or_B']
  #print 'INT_or_EXT after adding fourth edge = ', G.es["INT_or_EXT"]
  return G
#End
###################################

###################################
#########Remove a tadpole##########
###################################
#Begin
def remove_tad(G,M):
  l = L(M)
  #count_tadpole=0
  tad_list=[]
  for i in range (1,l):
    x=G.get_edgelist()[i]
    #print i, x
    if x[0]==x[1]:
      tad_list.append(x[0])
  count_tadpole=len(tad_list)
  print 'count_tadpole = ', count_tadpole
  print 'tadpol_list = ', tad_list
  if count_tadpole>0:
    tad_index = randint(0,count_tadpole-1)
    print 'tad_index = ', tad_index
    tad_node = tad_list[tad_index]
    print 'tad_node = ', tad_node
    tad_neigh = G.neighbors(tad_node)
    print tad_neigh
    G.delete_edges(G.get_eid(tad_node,tad_node))
    tad_neigh_in = G.neighbors(tad_node,'IN')
    tad_neigh_out = G.neighbors(tad_node,'OUT')
    print tad_neigh_in
    print tad_neigh_out
    print G
    if len(tad_neigh_out)!=0:
      G.delete_edges(G.get_eid(tad_node,tad_neigh_out[0]))
      V=tad_neigh_out[0]
    else:
      G.delete_edges(G.get_eid(tad_neigh_in[0],tad_node))
      V=tad_neigh_in[0]
    print 'tad_leg = ', V 
    V_con=[]
    V_neigh_in = G.neighbors(V,'IN')
    V_neigh_out = G.neighbors(V,'OUT')
    for i in range (0,len(V_neigh_out)):
      G.delete_edges(G.get_eid(V,V_neigh_out[i]))
      V_con.append(V_neigh_out[i])
    for i in range (0,len(V_neigh_in)):
      G.delete_edges(G.get_eid(V_neigh_in[i],V))
      V_con.append(V_neigh_in[i])
    print V_con
    print G
    G.add_edges([(V_con[1],V_con[0])])  # Add a fermionic line
    G.es[G.get_eid(V_con[1],V_con[0])]['F_or_B'] = 1   # Update 'F_or_B' attribute of the new edge
    G.es[G.get_eid(V_con[1],V_con[0])]['INT_or_EXT'] = 1   # Update 'INT_or_EXT' attribute of the new edge
    print G
    G.delete_vertices(tad_node)
    if tad_node<V:
      G.delete_vertices(V-1)
    else:
      G.delete_vertices(V) 
    print G
#End
###################################


#m=1
#g=generate_g_1_con(m)

m=2
g=generate_g_2(m)
print g
add_tad(g,m)
print g
add_tad(g,m+1)
print g
add_tad(g,m+2)
print g
remove_tad(g,m+3)

  
