##################################
##########Import Libraries########
##################################
#Begin
from Label_Pol import *
from igraph import *
from random import *
import numpy
#End
###################################

###################################
###########Add vertices############
###################################
#Begin
def add_nodes(G):
  NV = G.vcount()
  #print 'NV = ', NV
  l = G.ecount()
  #print 'l = ', l
  l_f = 0
  for i in range (0,l):   # Sweep all the lines
    if G.es[i]['F_or_B'] == 1:   # If the line is fermionic
      l_f = l_f+1  
  #print 'l_f = ', l_f
  f_list = [None]*l_f   # To store the fermionic lines' indices
  index_f_list = 0   # The index of f_list
  for i in range (0,l):   # Sweep all the lines
    if G.es[i]['F_or_B'] == 1:   # If the line is fermionic
      f_list[index_f_list] = i   # Store the fermionic line index
      index_f_list = index_f_list + 1 
  #print 'f_list = ', f_list  
  f1 = randint(0,l_f-1)   # Choose one of the fermionic lines randomly
  f2 = randint(0,l_f-1)   # Choose another fermionic line randomly
  #print 'f1 = ', f1
  #print 'f2 = ', f2
  V1 = G.get_edgelist()[f_list[f1]][0]  # The first node connected to f1
  #print 'V1 = ', V1
  V2 = G.get_edgelist()[f_list[f1]][1]  # The second node connected to f1
  #print 'V2 = ', V2
  V3 = G.get_edgelist()[f_list[f2]][0]  # The first node connected to f1
  #print 'V3 = ', V3
  V4 = G.get_edgelist()[f_list[f2]][1]  # The second node connected to f1
  #print 'V4 = ', V4
  #print 'Initial diagram = ', G
  #print 'Initial name1 = ', G.vs['name1']
  #print 'Initial visited = ', G.vs['visited']
  #print 'Initial F_or_B = ', G.es['F_or_B']
  #print 'Initial INT_or_EXT = ', G.es["INT_or_EXT"]
  if f1==f2:
    G.add_vertices(2)  # Add a node
    G.vs[NV]['name1'] = str(NV)  # Update 'name1' attribute of the new node
    G.vs[NV]["visited"] = 0       # Update 'visited' attribute of the new node
    G.vs[NV+1]['name1'] = str(NV+1)  # Update 'name1' attribute of the new node
    G.vs[NV+1]["visited"] = 0       # Update 'visited' attribute of the new node
    #print 'Diagram after one node = ', G
    #print 'name1 after one node = ', G.vs['name1']
    #print 'visited after one node = ', G.vs['visited']
    #print 'F_or_B after one node = ', G.es['F_or_B']
    #print 'INT_or_EXT after one node = ', G.es["INT_or_EXT"]
    G.delete_edges(f_list[f1])  # Delete f1
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
    G.add_edges([(NV,NV+1)])     # Add second edge between new node and V2
    G.es[G.get_eid(NV,NV+1)]['F_or_B'] = 1   # Update 'F_or_B' attribute of the new edge
    G.es[G.get_eid(NV,NV+1)]['INT_or_EXT'] = 1   # Update 'INT_or_EXT' attribute of the new edge
    #print 'Diagram after adding second edge = ', G
    #print 'name1 after adding second edge = ', G.vs['name1']
    #print 'visited after adding second edge = ', G.vs['visited']
    #print 'F_or_B after adding second edge = ', G.es['F_or_B']
    #print 'INT_or_EXT after adding second edge = ', G.es["INT_or_EXT"]
    G.add_edges([(NV+1,V2)])     # Add second edge between new node and V2
    G.es[G.get_eid(NV+1,V2)]['F_or_B'] = 1   # Update 'F_or_B' attribute of the new edge
    G.es[G.get_eid(NV+1,V2)]['INT_or_EXT'] = 1   # Update 'INT_or_EXT' attribute of the new edge
    #print 'Diagram after adding second edge = ', G
    #print 'name1 after adding second edge = ', G.vs['name1']
    #print 'visited after adding second edge = ', G.vs['visited']
    #print 'F_or_B after adding second edge = ', G.es['F_or_B']
    #print 'INT_or_EXT after adding second edge = ', G.es["INT_or_EXT"]
  if f1!=f2:
    G.add_vertices(2)  # Add two nodes
    G.vs[NV]['name1'] = str(NV)  # Update 'name1' attribute of the new node
    G.vs[NV]["visited"] = 0       # Update 'visited' attribute of the new node
    G.vs[NV+1]['name1'] = str(NV+1)  # Update 'name1' attribute of the new node
    G.vs[NV+1]["visited"] = 0       # Update 'visited' attribute of the new node
    #print 'Diagram after one node = ', G
    #print 'name1 after one node = ', G.vs['name1']
    #print 'visited after one node = ', G.vs['visited']
    #print 'F_or_B after one node = ', G.es['F_or_B']
    #print 'INT_or_EXT after one node = ', G.es["INT_or_EXT"]
    G.delete_edges(f_list[f1])  # Delete f1
    #print '*****************8'
    #print 'Diagram after deleting one edge = ', G
    #print 'name1 after deleting one edge = ', G.vs['name1']
    #print 'visited after deleting one edge = ', G.vs['visited']
    #print 'F_or_B after deleting one edge = ', G.es['F_or_B']
    #print 'INT_or_EXT after deleting one edge = ', G.es["INT_or_EXT"]
    if f_list[f1]<f_list[f2]:
      G.delete_edges(f_list[f2]-1)  # Delete f2
    else:
      G.delete_edges(f_list[f2])  # Delete f2
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
    
    G.add_edges([(V3,NV+1)])     # Add first edge between V1 and new node
    G.es[G.get_eid(V3,NV+1)]['F_or_B'] = 1   # Update 'F_or_B' attribute of the new edge
    G.es[G.get_eid(V3,NV+1)]['INT_or_EXT'] = 1   # Update 'INT_or_EXT' attribute of the new edge
    #print 'Diagram after adding first edge = ', G
    #print 'name1 after adding first edge = ', G.vs['name1']
    #print 'visited after adding first edge = ', G.vs['visited']
    #print 'F_or_B after adding first edge = ', G.es['F_or_B']
    #print 'INT_or_EXT after adding first edge = ', G.es["INT_or_EXT"]
    G.add_edges([(NV+1,V4)])     # Add second edge between new node and V2
    G.es[G.get_eid(NV+1,V4)]['F_or_B'] = 1   # Update 'F_or_B' attribute of the new edge
    G.es[G.get_eid(NV+1,V4)]['INT_or_EXT'] = 1   # Update 'INT_or_EXT' attribute of the new edge
    #print 'Diagram after adding second edge = ', G
    #print 'name1 after adding second edge = ', G.vs['name1']
    #print 'visited after adding second edge = ', G.vs['visited']
    #print 'F_or_B after adding second edge = ', G.es['F_or_B']
    #print 'INT_or_EXT after adding second edge = ', G.es["INT_or_EXT"]   
  return G, NV   
#End
################################### 

###################################
########Add interaction line#######
###################################
#Begin
def add_int_line(G):
  out1 = add_nodes(G)
  first_node = out1[1]
  second_node = first_node+1
  G_imt = out1[0]
  #print 'G_imt = ', G_imt
  #print 'First added node = ', first_node
  #print 'Second added node = ', second_node
  #print 'G_final_before_int_line = ', G_imt
  
  G.add_edges([(first_node,second_node)])     # Add an interaction line between two new nodes
  G.es[G.get_eid(first_node,second_node)]['F_or_B'] = 0   # Update 'F_or_B' attribute of the new edge
  G.es[G.get_eid(first_node,second_node)]['INT_or_EXT'] = 1   # Update 'INT_or_EXT' attribute of the new edge
  #print 'Final diagram = ', G
  #print 'name1 = ', G.vs['name1']
  #print 'visited = ', G.vs['visited']
  #print 'F_or_B = ', G.es['F_or_B']
  #print 'INT_or_EXT = ', G.es["INT_or_EXT"]
  return G
#End
###################################

###################################
######Remove Interaction line######
###################################
#Begin
def remove_int_line(G):
  #print 'Initial diagram = ', G
  #print 'Initial name1 = ', G.vs['name1']
  #print 'Initial visited = ', G.vs['visited']
  #print 'Initial F_or_B = ', G.es['F_or_B']
  #print 'Initial INT_or_EXT = ', G.es["INT_or_EXT"]
  
  NV = G.vcount()   # Total number od nodes
  #print 'NV = ', NV  
  l = G.ecount()   # Total number of lines
  #print 'l = ', l
  l_b = 0  # To store total number of interaction lines
  legs = bubble_finder(G)[1]
  #print '******#####*****'
  #print legs
  for i in range (0,l):   # Sweep all the lines
    if G.es[i]['F_or_B'] == 0 and G.es[i]['INT_or_EXT']==1:   # If the line is fermionic
      V1 = G.get_edgelist()[i][0]  # The first node connected to b
      V2 = G.get_edgelist()[i][1]  # The second node connected to b
      flag = True
      #for j in range(0,len(legs)):
        #if V1==legs[j] or V2==legs[j]:
          #flag = False
      if flag:
        l_b = l_b + 1
  if l_b>0:
    #print 'l_b = ', l_b
    b_list = [None]*l_b   # To store the bosonic lines' indices
    index_b_list = 0   # The index of b_list
    for i in range (0,l):   # Sweep all the lines
      if G.es[i]['F_or_B'] == 0 and G.es[i]['INT_or_EXT']==1:   # If the line is internal and bosonic
        V1 = G.get_edgelist()[i][0]  # The first node connected to b
        V2 = G.get_edgelist()[i][1]  # The second node connected to b
        flag = True
        #for j in range(0,len(legs)):
          #if V1==legs[j] or V2==legs[j]:
            #flag = False
        if flag:
          b_list[index_b_list] = i   # Store the bosonic line index
          index_b_list = index_b_list + 1 
    #print 'b_list = ', b_list  
    b = randint(0,l_b-1)   # Choose one of the bosonic lines randomly
    #print 'b = ', b
    V1 = G.get_edgelist()[b_list[b]][0]  # The first node connected to b
    #print 'V1 = ', V1
    V2 = G.get_edgelist()[b_list[b]][1]  # The second node connected to b
    #print 'V2 = ', V2
    
    G.delete_edges(b_list[b])  # Remove the interaction line
    
    #print 'Diagram after deleting one b-line = ', G
    #print 'name1 after deleting one b-line = ', G.vs['name1']
    #print 'visited after deleting one b-line = ', G.vs['visited']
    #print 'F_or_B after deleting one b-line = ', G.es['F_or_B']
    #print 'INT_or_EXT after deleting one b-line = ', G.es["INT_or_EXT"]
    
    n1_IN = G.neighbors(V1, mode="IN")   # Find ingoing neughbors of V1
    #print 'In neighbors of V1 =', n1_IN
    #print 'V1 = ', V1
    
    
    for i in range (0,len(n1_IN)):
      G.delete_edges(G.get_eid(n1_IN[i],V1))  # Remove adjacent ingoing edges of V1
    #print 'Diagram after deleting f-lines connected to V1 = ', G
    #print 'name1 after deleting f-lines connected to V1 = ', G.vs['name1']
    #print 'visited after deleting f-lines connected to V1 = ', G.vs['visited']
    #print 'F_or_B after deleting f-lines connected to V1 = ', G.es['F_or_B']
    #print 'INT_or_EXT after deleting f-lines connected to V1 = ', G.es["INT_or_EXT"]

    n1_OUT = G.neighbors(V1, mode="OUT")
    #print 'OUT neighbors of V1 =', n1_OUT
    for i in range (0,len(n1_OUT)):
      G.delete_edges(G.get_eid(V1,n1_OUT[i]))  # Remove adjacent ingoing edges of V1
    #print 'Diagram after deleting f-lines connected to V1 = ', G
    #print 'name1 after deleting f-lines connected to V1 = ', G.vs['name1']
    #print 'visited after deleting f-lines connected to V1 = ', G.vs['visited']
    #print 'F_or_B after deleting f-lines connected to V1 = ', G.es['F_or_B']
    #print 'INT_or_EXT after deleting f-lines connected to V1 = ', G.es["INT_or_EXT"]
    v = [None]*2
    j = 0
    for i in range (0,len(n1_IN)):
      v[j] = n1_IN[i]
      j = j+1
    for i in range (0,len(n1_OUT)):
      v[j] = n1_OUT[i]
      j = j+1
    G.add_edges([(v[0],v[1])])  # Add a fermionic line
    G.es[G.get_eid(v[0],v[1])]['F_or_B'] = 1   # Update 'F_or_B' attribute of the new edge
    G.es[G.get_eid(v[0],v[1])]['INT_or_EXT'] = 1   # Update 'INT_or_EXT' attribute of the new edge
    #print G


    n2_IN = G.neighbors(V2, mode="IN")   # Find ingoing neughbors of V1
    #print 'In neighbors of V2 =', n1_IN
    #print 'V2 = ', V2
    
    
    for i in range (0,len(n2_IN)):
      G.delete_edges(G.get_eid(n2_IN[i],V2))  # Remove adjacent ingoing edges of V1
    #print 'Diagram after deleting f-lines connected to V2 = ', G
    #print 'name1 after deleting f-lines connected to V2 = ', G.vs['name1']
    #print 'visited after deleting f-lines connected to V2 = ', G.vs['visited']
    #print 'F_or_B after deleting f-lines connected to V2 = ', G.es['F_or_B']
    #print 'INT_or_EXT after deleting f-lines connected to V2 = ', G.es["INT_or_EXT"]

    n2_OUT = G.neighbors(V2, mode="OUT")
    #print 'OUT neighbors of V2 =', n2_OUT
    for i in range (0,len(n2_OUT)):
      G.delete_edges(G.get_eid(V2,n2_OUT[i]))  # Remove adjacent ingoing edges of V1
    #print 'Diagram after deleting f-lines connected to V2 = ', G
    #print 'name1 after deleting f-lines connected to V2 = ', G.vs['name1']
    #print 'visited after deleting f-lines connected to V2 = ', G.vs['visited']
    #print 'F_or_B after deleting f-lines connected to V2 = ', G.es['F_or_B']
    #print 'INT_or_EXT after deleting f-lines connected to V2 = ', G.es["INT_or_EXT"]
    v = [None]*2
    j = 0
    for i in range (0,len(n2_IN)):
      v[j] = n2_IN[i]
      j = j+1
    for i in range (0,len(n2_OUT)):
      v[j] = n2_OUT[i]
      j = j+1
    G.add_edges([(v[0],v[1])])  # Add a fermionic line
    G.es[G.get_eid(v[0],v[1])]['F_or_B'] = 1   # Update 'F_or_B' attribute of the new edge
    G.es[G.get_eid(v[0],v[1])]['INT_or_EXT'] = 1   # Update 'INT_or_EXT' attribute of the new edge
    #print G
    G.delete_vertices(V1)
    if V1<V2:
      G.delete_vertices(V2-1)
    else:
      G.delete_vertices(V2) 
    for i in range (0,G.vcount()):
      G.vs[i]['name1'] = i
    #print 'Final Diagram = ', G
    #print 'Final name1 = ', G.vs['name1']
    #print 'Final visited = ', G.vs['visited']
    #print 'Final F_or_B = ', G.es['F_or_B']
    #print 'Final INT_or_EXT = ', G.es["INT_or_EXT"]  
  return G
#End
###################################

###################################
#########Bubble generation#########
###################################
#Begin
def add_bubble(G):
  #print 'Initial diagram = ', G
  #print 'Initial name1 = ', G.vs['name1']
  #print 'Initial visited = ', G.vs['visited']
  #print 'Initial F_or_B = ', G.es['F_or_B']
  #print 'Initial INT_or_EXT = ', G.es["INT_or_EXT"]
  
  NV = G.vcount()   # Total number od nodes
  #print 'NV = ', NV  
  l = G.ecount()   # Total number of lines
  #print 'l = ', l
  l_b = 0  # To store total number of internal interaction lines
  for i in range (0,l):   # Sweep all the lines
    if G.es[i]['F_or_B'] == 0: #and G.es[i]['INT_or_EXT']==1:   # If the line is fermionic
      l_b = l_b+1  
  if l_b>=0:  # if True
    #print 'l_b = ', l_b
    b_list = [None]*(l_b)   # To store the bosonic lines' indices
    index_b_list = 0   # The index of b_list
    for i in range (0,l):   # Sweep all the lines
      if G.es[i]['F_or_B'] == 0: #and G.es[i]['INT_or_EXT']==1:   # If the line is internal and bosonic
        b_list[index_b_list] = i   # Store the bosonic line index
        index_b_list = index_b_list + 1 
    #print 'b_list = ', b_list  
    #b = randint(0,l_b-1)   # Choose one of the internal bosonic lines randomly
    b = randint(0,l_b-1)   # Choose one of the internal bosonic lines randomly
    #print 'b = ', b
    
    V1 = G.get_edgelist()[b_list[b]][0]  # The first node connected to b
    #print 'V1 = ', V1
    V2 = G.get_edgelist()[b_list[b]][1]  # The second node connected to b
    #print 'V2 = ', V2
    
    G.delete_edges(b_list[b])  # Remove the interaction line
    
    #print 'Diagram after deleting one b-line = ', G
    #print 'name1 after deleting one b-line = ', G.vs['name1']
    #print 'visited after deleting one b-line = ', G.vs['visited']
    #print 'F_or_B after deleting one b-line = ', G.es['F_or_B']
    #print 'INT_or_EXT after deleting one b-line = ', G.es["INT_or_EXT"]
    G.add_vertices(2)  # Add a node
    G.vs[NV]['name1'] = str(NV)  # Update 'name1' attribute of the new node
    G.vs[NV]["visited"] = 0       # Update 'visited' attribute of the new node
    G.vs[NV+1]['name1'] = str(NV+1)  # Update 'name1' attribute of the new node
    G.vs[NV+1]["visited"] = 0       # Update 'visited' attribute of the new node
    #print 'Diagram after one node = ', G
    #print 'name1 after one node = ', G.vs['name1']
    #print 'visited after one node = ', G.vs['visited']
    #print 'F_or_B after one node = ', G.es['F_or_B']
    #print 'INT_or_EXT after one node = ', G.es["INT_or_EXT"]
    
    G.add_edges([(V1,NV)])     # Add an interaction line between two new nodes
    G.es[G.get_eid(V1,NV)]['F_or_B'] = 0   # Update 'F_or_B' attribute of the new edge
    #G.es[G.get_eid(V1,NV)]['INT_or_EXT'] = 1   # Update 'INT_or_EXT' attribute of the new edge

    G.add_edges([(NV+1,V2)])     # Add an interaction line between two new nodes
    G.es[G.get_eid(NV+1,V2)]['F_or_B'] = 0   # Update 'F_or_B' attribute of the new edge
    #G.es[G.get_eid(NV+1,V2)]['INT_or_EXT'] = 1   # Update 'INT_or_EXT' attribute of the new edge

    G.add_edges([(NV,NV+1)])     # Add first fermionic edge
    G.es[G.get_eid(NV,NV+1)]['F_or_B'] = 1   # Update 'F_or_B' attribute of the new edge
    G.es[G.get_eid(NV,NV+1)]['INT_or_EXT'] = 1   # Update 'INT_or_EXT' attribute of the new edge

    G.add_edges([(NV+1,NV)])     # Add first fermionic edge
    G.es[G.get_eid(NV+1,NV)]['F_or_B'] = 1   # Update 'F_or_B' attribute of the new edge
    G.es[G.get_eid(NV+1,NV)]['INT_or_EXT'] = 1   # Update 'INT_or_EXT' attribute of the new edge

    if (len(G.adjacent(V1, mode=IN)) + len(G.adjacent(V1, mode=OUT)) == 3):
      G.es[G.get_eid(V1,NV)]['INT_or_EXT'] = 1 
    else:
      G.es[G.get_eid(V1,NV)]['INT_or_EXT'] = 0 
    if (len(G.adjacent(V2, mode=IN)) + len(G.adjacent(V2, mode=OUT)) == 3):
      G.es[G.get_eid(NV+1,V2)]['INT_or_EXT'] = 1 
    else:
      G.es[G.get_eid(NV+1,V2)]['INT_or_EXT'] = 0 

    #print 'Diagram after adding second edge = ', G
    #print 'name1 after adding second edge = ', G.vs['name1']
    #print 'visited after adding second edge = ', G.vs['visited']
    #print 'F_or_B after adding second edge = ', G.es['F_or_B']
    #print 'INT_or_EXT after adding second edge = ', G.es["INT_or_EXT"] 
  return G  
#End
###################################

###################################
#########Bubble Finder############# 
###################################
#Begin
def bubble_finder(G):
  #print 'Initial diagram = ', G
  #print 'Initial name1 = ', G.vs['name1']
  #print 'Initial visited = ', G.vs['visited']
  #print 'Initial F_or_B = ', G.es['F_or_B']
  #print 'Initial INT_or_EXT = ', G.es["INT_or_EXT"]
  NV = G.vcount()   # Total number od nodes
  #print 'NV = ', NV  
  l = G.ecount()   # Total number of lines
  #print 'l = ', l
  bubble_list = []
  nodes_list = []
  legs = []
  
  for i in range(0,NV):
    #print i
    counter = 0
    out_nodes = G.neighbors(i,mode=OUT)
    in_nodes = G.neighbors(i,mode=IN)
    nodes = [None]*(len(out_nodes) + len(in_nodes))
    for j in range(0,len(out_nodes)):
      nodes[counter]=out_nodes[j]
      counter = counter + 1
    for j in range(0,len(in_nodes)):
      nodes[counter]=in_nodes[j]
      counter = counter + 1
    #print 'len = ', len(nodes)
    #print 'node = ', i
    #print 'out_nodes = ', out_nodes
    #print 'in_nodes = ', in_nodes
    for j in range(0,len(nodes)):
      for k in range(j+1,len(nodes)): 
        #print 'j =',j
        #print 'k =', k
        a = nodes[j]
        b = nodes[k]
        
        if a == b:
          
          if a in in_nodes:
            if a in out_nodes:
              if G.es[G.get_eid(a,i)]['F_or_B'] == 1 and G.es[G.get_eid(i,a)]['F_or_B'] == 1 and a!=i:
                bubble_list.append(nodes)
                nodes_list.append(a)
                nodes_list.append(i)
                for r in range (0,len(nodes)):
                  if nodes[r]!=a:
                    legs.append(nodes[r])
                for s in range(0,len(G.neighbors(a,mode=OUT))):
                  if G.es[G.get_eid(a,G.neighbors(a,mode=OUT)[s])]['F_or_B'] == 0:
                    legs.append(G.neighbors(a,mode=OUT)[s])
                for s in range(0,len(G.neighbors(a,mode=IN))):
                  if G.es[G.get_eid(G.neighbors(a,mode=IN)[s],a)]['F_or_B'] == 0:
                    legs.append(G.neighbors(a,mode=IN)[s])  
          
          if a in in_nodes:
            if a not in out_nodes:
              if len(out_nodes)!=0:
                c = out_nodes[0] 
              else:
                for t in range(0,len(in_nodes)):
                  if in_nodes[t]!=a:
                    c=in_nodes[t] 
              if (c in in_nodes) and (G.es[G.get_eid(c,i)]['F_or_B'] == 0) and (a!=i):
                bubble_list.append(nodes)
                nodes_list.append(a)
                nodes_list.append(i)
                for r in range (0,len(nodes)):
                  if nodes[r]!=a:
                    legs.append(nodes[r])
                for s in range(0,len(G.neighbors(a,mode=OUT))):
                  if G.es[G.get_eid(a,G.neighbors(a,mode=OUT)[s])]['F_or_B'] == 0:
                    legs.append(G.neighbors(a,mode=OUT)[s])
                for s in range(0,len(G.neighbors(a,mode=IN))):
                  if G.es[G.get_eid(G.neighbors(a,mode=IN)[s],a)]['F_or_B'] == 0:
                    legs.append(G.neighbors(a,mode=IN)[s])  
              if (c in out_nodes) and (G.es[G.get_eid(i,c)]['F_or_B'] == 0) and (a!=i):
                bubble_list.append(nodes)
                nodes_list.append(a)
                nodes_list.append(i)
                for r in range (0,len(nodes)):
                  if nodes[r]!=a:
                    legs.append(nodes[r])
                for s in range(0,len(G.neighbors(a,mode=OUT))):
                  if G.es[G.get_eid(a,G.neighbors(a,mode=OUT)[s])]['F_or_B'] == 0:
                    legs.append(G.neighbors(a,mode=OUT)[s])
                for s in range(0,len(G.neighbors(a,mode=IN))):
                  if G.es[G.get_eid(G.neighbors(a,mode=IN)[s],a)]['F_or_B'] == 0:
                    legs.append(G.neighbors(a,mode=IN)[s])    

          if a not in in_nodes:
            if a in out_nodes:
              if len(in_nodes)!=0:
                c = in_nodes[0] 
              else:
                for t in range(0,len(out_nodes)):
                  if out_nodes[t]!=a:
                    c=out_nodes[t] 
              if (c in in_nodes) and (G.es[G.get_eid(c,i)]['F_or_B'] == 0) and (a!=i):
                bubble_list.append(nodes)
                nodes_list.append(a)
                nodes_list.append(i)
                for r in range (0,len(nodes)):
                  if nodes[r]!=a:
                    legs.append(nodes[r])
                for s in range(0,len(G.neighbors(a,mode=OUT))):
                  if G.es[G.get_eid(a,G.neighbors(a,mode=OUT)[s])]['F_or_B'] == 0:
                    legs.append(G.neighbors(a,mode=OUT)[s])
                for s in range(0,len(G.neighbors(a,mode=IN))):
                  if G.es[G.get_eid(G.neighbors(a,mode=IN)[s],a)]['F_or_B'] == 0:
                    legs.append(G.neighbors(a,mode=IN)[s])  
              if (c in out_nodes) and (G.es[G.get_eid(i,c)]['F_or_B'] == 0) and (a!=i):
                bubble_list.append(nodes)
                nodes_list.append(a)
                nodes_list.append(i)
                for r in range (0,len(nodes)):
                  if nodes[r]!=a:
                    legs.append(nodes[r])
                for s in range(0,len(G.neighbors(a,mode=OUT))):
                  if G.es[G.get_eid(a,G.neighbors(a,mode=OUT)[s])]['F_or_B'] == 0:
                    legs.append(G.neighbors(a,mode=OUT)[s])
                for s in range(0,len(G.neighbors(a,mode=IN))):
                  if G.es[G.get_eid(G.neighbors(a,mode=IN)[s],a)]['F_or_B'] == 0:
                    legs.append(G.neighbors(a,mode=IN)[s])  
  
    #for j in range (0, len(nodes)):
  #print bubble_list  
  #print nodes_list
  #print legs
  for i in range(0,len(legs)/2):
    for j in range(i+1,len(legs)/2):
      if (legs[2*i]==legs[2*j] and legs[2*i+1]==legs[2*j+1]) or (legs[2*i]==legs[2*j+1] and legs[2*i+1]== legs[2*j]):
        legs[2*j]=None
        legs[2*j+1]=None
        nodes_list[2*j]=None
        nodes_list[2*j+1]=None
  new_legs=[]
  new_nodes_list=[]
  for i in range(0,len(legs)):
    if legs[i]!=None:
      new_legs.append(legs[i])
      new_nodes_list.append(nodes_list[i])
  return new_nodes_list, new_legs
#End
###################################

###################################
###########Remove bubble###########
###################################
#Begin
def remove_bubble(G):
  bubbles = bubble_finder(G)  # Find bubbles
  bubble_nodes = bubbles[0]   # Bubbles' nodes
  bubble_legs = bubbles[1]    # Bubbles' legs
  n_bubble = len(bubble_nodes)/2
  #print bubble_nodes
  #print bubble_legs
  #print 'number of bubbles = ', n_bubble
  if n_bubble!=0:
    bubble_choice = randint(0, n_bubble-1) 
    #print 'bubble_choice = ', bubble_choice
    Vb_1 = bubble_nodes[2*bubble_choice]
    Vb_2 = bubble_nodes[2*bubble_choice+1]
    Vl_1 = bubble_legs[2*bubble_choice]
    Vl_2 = bubble_legs[2*bubble_choice+1]
    #print 'bubble nodes = ', Vb_1, Vb_2
    #print 'bubble legs = ', Vl_1, Vl_2
    G.add_edges([(Vl_1,Vl_2)])     # Add an interaction line between two legs
    G.es[G.get_eid(Vl_1,Vl_2)]['F_or_B'] = 0   # Update 'F_or_B' attribute of the new edge
    G.es[G.get_eid(Vl_1,Vl_2)]['INT_or_EXT'] = 1   # Update 'INT_or_EXT' attribute of the new edge
    G.delete_vertices(Vb_1)
    if Vb_1<Vb_2:
      G.delete_vertices(Vb_2-1) 
    else:
      G.delete_vertices(Vb_2)
    for i in range (0,G.vcount()):
      G.vs[i]['name1']=str(i)
  return G
#End
###################################


