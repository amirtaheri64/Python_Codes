#######################################
##########Import Libraries#############
#######################################
#Begin
from igraph import *
from random import *
from Hubbard_like_diagram import find_prc_line
import numpy
#End
#######################################

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
  flag = False
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
##########Connected_loops##############
#######################################
#Begin
def connected_loops(G,V):
  loop = []
  loop.append(loop_finder(G,V))

  #print 'loop = ', loop
  connections_out = []
  loop_out = []
  connections_in = []
  loop_in = []
  for i in loop[0]:
    x_out = G.neighbors(i,mode="OUT")
    for j in range (0,len(x_out)):
      if x_out[j] not in loop[0]:
        connections_out.append(i) 
        if G.vs[x_out[j]]["visited"] == 0:
          loop_out.append(loop_finder(G,x_out[j])) 
      
    x_in = G.neighbors(i,mode="IN")
    for j in range (0,len(x_in)):
      if x_in[j] not in loop[0]:
        connections_in.append(i) 
        if G.vs[x_in[j]]["visited"] == 0:
          loop_in.append(loop_finder(G,x_in[j])) 
    for j in range(0,len(loop_out)):
      G.vs[loop_out[j]]["visited"] = 0  
    for j in range(0,len(loop_in)):
      G.vs[loop_in[j]]["visited"] = 0 
  
  return connections_out,loop_out,connections_in,loop_in,loop 
#End
#######################################

#######################################
########Is G a connected graph#########
#######################################
#Begin
def is_connected(G):
  NV = G.vcount()   # Total number of nodes
  #print 'NV = ', NV  
  l = G.ecount()   # Total number of lines
  #print 'l = ', l
  #V=-1
  #while V==-1:
    #V = randint(0,NV-1)
    #print V
    #if len(G.neighbors(V,mode="OUT"))+len(G.neighbors(V,mode="IN"))==1:
      #V=-1
    #print 'node = ', V
  #print 'prc_line = ', find_prc_line(G)
  V = find_prc_line(G)[1]
  #print "V = ",V
  a = connected_loops(G,V)
  CONNECTIONS_OUT = a[0]
  LOOP_OUT = a[1]
  CONNECTIONS_IN = a[2]
  LOOP_IN = a[3]
  INITIAL_LOOP = a[4]
  #print 'connections_in = ', CONNECTIONS_IN
  #print 'loops_in = ', LOOP_IN
  #print 'connections_out = ', CONNECTIONS_OUT
  #print 'loops_out = ', LOOP_OUT

  LOOPS = INITIAL_LOOP+LOOP_IN + LOOP_OUT
  #print 'loops = ', LOOPS
  new = LOOPS
  if len(new[0])==0:
    return False
  s=1
  while s!=0:# and len(LOOP_IN)!=0 and len(LOOP_OUT)!=0:
    LOOPS = new
    #print 'new = ', new
    s=0
    for i in range (0,len(LOOPS)):
      #for j in range (0,len(LOOPS[i])):
        V = LOOPS[i][0]
        #print 'V = ', V
        a = connected_loops(G,V)
        CONNECTIONS_OUT = a[0]
        LOOP_OUT = a[1]
        CONNECTIONS_IN = a[2]
        LOOP_IN = a[3]
        s = s + len(CONNECTIONS_OUT) + len(CONNECTIONS_IN) 
        #print 'connections_in = ', CONNECTIONS_IN
        #print 'loops_in = ', LOOP_IN
        #print 'connections_out = ', CONNECTIONS_OUT
        #print 'loops_out = ', LOOP_OUT
        new = new + LOOP_IN + LOOP_OUT
      
  #print new  
  connected_nodes =[]
  for i in range(0,len(new)):
    for j in range(0,len(new[i])):
      if new[i][j] not in connected_nodes:
        connected_nodes.append(new[i][j])
  #print connected_nodes
  #print len(connected_nodes)
  #print NV

  if len(connected_nodes)==NV-2:
    #print 'This is a connected diagram!'
    val = True
  else:
    #print 'This is a disconnected diagram!'
    val = False
  for i in range(0,NV):
    G.vs[i]["visited"] = 0 
    if len(G.neighbors(i))==1:
      G.vs[i]["visited"] = 1   
  return val
#End
#######################################

  

