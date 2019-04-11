'''
Last update: March-12-2019
'''

##################################
##########Import Libraries########
##################################
#Begin
from igraph import *
from random import *
from Label_self import *
import numpy
#End
###################################

'''
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
  G.es["name2"] = ["0-F", "1-B", "2-F", "3-F", "4-B", "5-F", "6-F"]
  G.es["F_or_B"] = [1,0,1,1,0,1,1]
  G.es["INT_or_EXT"] = [0,1,1,1,1,1,0]
  G.es["Label"] = [[None]*(M+1+1)]*l  # Define "None" labels
   
  #End
  #####################################
  
  return G

#End
######################################

########################################
##Define a 3rd order self-energy graph##
########################################
#Begin
def generate_g_3(M):  # M: order of diagram
  n_v = N_V(M)
  l = L(M)
  G = Graph(directed=True)
  G.add_vertices(n_v) 

  G.add_edges([(0,2)])
  G.add_edges([(2,3)])
  G.add_edges([(3,4)])
  G.add_edges([(4,5)])
  G.add_edges([(5,6)])
  G.add_edges([(6,3)])
  G.add_edges([(4,6)])
  G.add_edges([(5,7)])
  G.add_edges([(2,7)])
  G.add_edges([(7,1)])

  #####################################
  ############Attributes###############
  #####################################
  #Begin
  G.vs["name1"] = [0, 1, 2, 3, 4, 5, 6, 7]
  G.vs["visited"] = [1,1,0,0,0,0,0,0]
  G.vs["Spin"] = [None]*n_v
  G.es["name2"] = ["0-F", "1-B", "2-F", "3-F", "4-F", "5-F", "6-B", "7-B", "8-F", "9-F"]
  G.es["F_or_B"] = [1,0,1,1,1,1,0,0,1,1]
  G.es["INT_or_EXT"] = [0,1,1,1,1,1,1,1,1,0]
  G.es["Label"] = [[None]*(M+1+1)]*l  # Define "None" labels
   
  #End
  #####################################
  
  return G

#End
######################################

########################################
##Define a 3rd order self-energy graph##
########################################
#Begin
def generate_g_3_tadpole(M):  # M: order of diagram
  n_v = N_V(M)
  l = L(M)
  G = Graph(directed=True)
  G.add_vertices(n_v) 

  G.add_edges([(0,2)])
  G.add_edges([(2,3)])
  G.add_edges([(3,4)])
  G.add_edges([(4,5)])
  G.add_edges([(5,3)])
  G.add_edges([(5,6)])
  G.add_edges([(2,6)])
  G.add_edges([(4,7)])
  G.add_edges([(7,7)])
  G.add_edges([(6,1)])

  #####################################
  ############Attributes###############
  #####################################
  #Begin
  G.vs["name1"] = [0, 1, 2, 3, 4, 5, 6, 7]
  G.vs["visited"] = [1,1,0,0,0,0,0,0]
  G.vs["Spin"] = [None]*n_v
  G.es["name2"] = ["0-F", "1-B", "2-F", "3-F", "4-F", "5-B", "6-F", "7-B", "8-F", "9-F"]
  G.es["F_or_B"] = [1,0,1,1,1,0,1,0,1,1]
  G.es["INT_or_EXT"] = [0,1,1,1,1,1,1,1,1,0]
  G.es["Label"] = [[None]*(M+1+1)]*l  # Define "None" labels
   
  #End
  #####################################
  
  return G

#End
######################################


########################################
##Define a 4th order self-energy graph##
########################################
#Begin
def generate_g_4(M):  # M: order of diagram
  n_v = N_V(M)
  l = L(M)
  G = Graph(directed=True)
  G.add_vertices(n_v) 

  G.add_edges([(0,2)])
  G.add_edges([(2,3)])
  G.add_edges([(3,4)])
  G.add_edges([(4,3)])
  G.add_edges([(4,5)])
  G.add_edges([(6,5)])
  G.add_edges([(7,6)])
  G.add_edges([(2,7)])
  G.add_edges([(7,8)])
  G.add_edges([(8,9)])
  G.add_edges([(9,8)])
  G.add_edges([(9,6)])
  G.add_edges([(5,1)])

  #####################################
  ############Attributes###############
  #####################################
  #Begin
  G.vs["name1"] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
  G.vs["visited"] = [1,1,0,0,0,0,0,0,0,0]
  G.vs["Spin"] = [None]*n_v
  G.es["name2"] = ["0-F", "1-B", "2-F", "3-F", "4-B", "5-F", "6-F", "7-F", "8_B", "9-F", "10-F", "11-B", "12-F"]
  G.es["F_or_B"] = [1,0,1,1,0,1,1,1,0,1,1,0,1]
  G.es["INT_or_EXT"] = [0,1,1,1,1,1,1,1,1,1,1,1,0]
  G.es["Label"] = [[None]*(M+1+1)]*l  # Define "None" labels
   
  #End
  #####################################
  
  return G

#End
######################################
'''
########################################
##Define a 2nd order self-energy graph##
########################################
#Begin
def generate_g_check(M):  # M: order of diagram
  n_v = N_V(M)
  l = L(M)
  G = Graph(directed=True)
  G.add_vertices(n_v) 

  G.add_edges([(0,4)])
  G.add_edges([(4,5)])
  G.add_edges([(4,5)])
  G.add_edges([(5,3)])
  G.add_edges([(2,3)])
  G.add_edges([(2,2)])
  G.add_edges([(3,1)])

  #####################################
  ############Attributes###############
  #####################################
  #Begin
  G.vs["name1"] = [0, 1, 2, 3, 4, 5]
  G.vs["visited"] = [1,1,0,0,0,0]
  G.vs["Spin"] = [None]*n_v
  G.es["name2"] = ["0-F", "1-F", "2-B", "3-F", "4-B", "5-F", "6-F"]
  G.es["F_or_B"] = [1,1,0,1,0,1,1]
  G.es["INT_or_EXT"] = [0,1,1,1,1,1,0]
  G.es["Label"] = [[None]*(M+1)]*l  # Define "None" labels
   
  #End
  #####################################
  
  return G

#End
######################################


##############################
#####Find principle line######
##############################
#Begin
def find_prc_line(G):
  order = 0
  for i in range (0,len(G.es["F_or_B"])):
    if G.es["F_or_B"][i] == 0:
      order=order+1 
  #print 'order = ', order
  node_in_out=[None]*2
  count_node_in_out=0
  for node in range (0,len(G.vs["visited"])):
    if G.vs["visited"][node] == 1:
      node_in_out[count_node_in_out] = node
      count_node_in_out=count_node_in_out+1
  #print 'in_out nodes = ', node_in_out 
  if len(G.neighbors(node_in_out[0], mode="OUT"))!=0:
    node_in = node_in_out[0]
    node_out = node_in_out[1]
  else:
    node_in = node_in_out[1]
    node_out = node_in_out[0]
  #print 'node_in = ', node_in
  #print 'node_out = ', node_out
  node = node_in
  #print 'node = ', node
  prc_nodes = [None]*N_V(order)
  #print 'prc_nodes = ', prc_nodes
  count_prc_nodes = 0
  while node!=node_out: 
    flag = True  
    prc_nodes[count_prc_nodes] = node
    count_prc_nodes = count_prc_nodes + 1 
    x = G.neighbors(node, mode="OUT")  # Out neighbors
    y=G.adjacent(node,OUT)
    #print 'out edges = ', y
    #print 'nodes_from_to = ', G.get_edgelist()[y[0]]
    #print '******', G.es[y[0]]["F_or_B"]
    for i in range (0,len(y)):
      if G.es[y[i]]["F_or_B"] == 1:  # If the outgoing line is fermionic
        #node = G.neighbors(node, mode="OUT")[i] # Go to the next node
        node = G.get_edgelist()[y[i]][1]
        y=G.adjacent(node,OUT)
        flag = False
        #print 'node = ', node
        #print flag
        break
    '''
    if flag:
      x = G.neighbors(node, mode="IN")  # Out neighbors
      for i in range (0,len(x)):
        if G.es[G.get_eid(node,x[i])]["F_or_B"] == 1:  # If the outgoing line is fermionic
          node = G.neighbors(node, mode="IN")[i] # Go to the next node
          print 'node = ', node
          break
     '''
  prc_nodes[count_prc_nodes] = node
  prc_nodes_fin = [None]*(count_prc_nodes+1)
  for i in range (0, count_prc_nodes+1):
    prc_nodes_fin[i] = prc_nodes[i]  
  #print 'prc_nodes = ', prc_nodes
  #print 'prc_nodes_fin = ', prc_nodes_fin
  #for i in range (0,count_prc_nodes):
    #print G.get_eid(prc_nodes_fin[i],prc_nodes_fin[i+1])
  #print order  
  return prc_nodes_fin

#######################################
#########loop finder procedure#########
#######################################
#Begin
def loop_finder(G,node):
  end_loop=0   # An auxuliay variable to 
  loop=[]
  #print 'node = ', node
  initial_node = node
  if G.vs[node]["visited"] == 0:
    loop.append(node)
  flag=False
  while not flag:
    
    G.vs[node]["visited"] = 1
    end_loop = 0
    flag = True
    y_out = G.adjacent(node, mode="OUT")
    y_in = G.adjacent(node, mode="IN") 
    x = G.neighbors(node,mode="OUT")  # Out neighbors
    #for i in range(0,len(y_out)):
      #print 'OUT_EDGE = ', G.get_edgelist()[y_out[i]]
    #for i in range(0,len(y_in)):
      #print 'IN_EDGE = ', G.get_edgelist()[y_in[i]]
    #print 'node = ', node
    for i in range (0,len(y_out)):
      x = G.get_edgelist()[y_out[i]][1]
      if G.es[y_out[i]]["F_or_B"] == 1 and G.vs[x]["visited"] == 0:  # If the outgoing line is fermionic
        node = x # Go to the next node
        #print 'new', node
        G.vs[node]["visited"] = 1
        loop.append(node)
        flag = False  
        break
    if flag:
      #print 'node = ', node 
      for j in range (0,len(y_in)):
        x = G.get_edgelist()[y_in[j]][0]
        if G.es[y_in[j]]["F_or_B"] == 1 and G.vs[x]["visited"] == 0:  # If the ingoing line is fermionic
          node = x  # Go to the next node
          #print 'new', node
          G.vs[node]["visited"] = 1
          loop.append(node)
          flag = False
         
          break    
    #print loop  
    #if flag==True:
      #end_loop=2
    
    #print 'Diagram after deleting one b-line = ', G
    #print 'name1 after deleting one b-line = ', G.vs['name1']
    #print 'visited after deleting one b-line = ', G.vs['visited']
    #print 'F_or_B after deleting one b-line = ', G.es['F_or_B']
    #print 'INT_or_EXT after deleting one b-line = ', G.es["INT_or_EXT"]  
  #print 'loop = ', loop 
   
  return loop
#End
#######################################

#######################################
##########Find all the loops###########
#######################################
#Begin
def all_loops(G,M):
  n_v = N_V(M) # number of vertices
  principle_line = find_prc_line(G) 

  #print 'principle line = ', principle_line

  all_nodes = principle_line

  loops=[]

  loops.append(principle_line)

  for i in range (0,n_v): # Sweep all the vertices
    if i not in all_nodes:
      x = loop_finder(G,i)
      loops.append(x)
      all_nodes=all_nodes+x
    if len(all_nodes) == n_v:
      i=n_v+1
  #print 'loops = ', loops
  #print 'all_nodes = ', all_nodes
  return loops
#End
####################################

def two_powers(n):
  val = 1
  for i in range(0,n):
    val = 2*val
  return val

def dec_to_bin(n):
  a=[]
  while n >1:
    a.append(n%2)
    n = n/2
  a.append(n)
  return a


####################################
###########Assign Spin##############
####################################
def Hubbard_diagram(G,M):
  flag = False
  loops=all_loops(G,M)
  #print 'loops = ', loops
  for i in range(0,len(loops[0])):
    G.vs[loops[0][i]]["Spin"] = 1
  #print G.vs["Spin"]
  for i in range (0,two_powers(len(loops)-1)):
    spins = dec_to_bin(i)
    for j in range(0, len(loops)-len(spins)-1):
      spins.append(0)
    #print spins
    
    for j in range (1,len(loops)):
      for k in range (0,len(loops[j])):
        G.vs[loops[j][k]]["Spin"] = spins[len(spins)-j] 
    #print G.vs["Spin"]
    check = 0   
    for j in range (0,len(G.es["F_or_B"])):
      if G.es["F_or_B"][j] == 0:
        b_line = G.get_edgelist()[j] 
        if G.vs["Spin"][b_line[0]]==G.vs["Spin"][b_line[1]]:
          check = check + 1
          j=len(G.es["F_or_B"])+1
    if check == 0:
      flag = True
      i = two_powers(len(loops)-1) + 1
      break
  return flag, G.vs["Spin"] 
#End
######################################

'''
g=generate_g_check(2)
print 'nodes tags =', g.vs["name1"]
print 'visited = ', g.vs["visited"] 
print 'spin = ', g.vs["Spin"]  
print 'F_or_B = ',  g.es["F_or_B"]
print 'INT or EXT = ', g.es["INT_or_EXT"] 
visual_style = {} 
color_dict = {"m": "black", "f": "white"}
visual_style["vertex_size"] = 20
visual_style["vertex_label"] = g.vs["name1"]
visual_style["edge_label"] = g.es["F_or_B"]
plot(g, **visual_style)  # Plot g

y=g.adjacent(4,OUT)
print 'y = ', y
print g.get_edgelist()[y[0]]
node = g.get_edgelist()[y[0]][1]
print node
print g.es[y[1]]["F_or_B"]

print find_prc_line(g)
'''

