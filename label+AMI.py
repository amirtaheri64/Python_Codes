'''
label+AMI: Get a diagram, label it and generate AMI arrays

Last updates: April-05-2019

'''
###################################
##########Import libraries#########
###################################
#Begin
from Label_self import *
from Hubbard_like_diagram import *
from is_connected import *
from Symbolic_multi_AMI_new import *
from add_remove import *
from AMI_zero import *
#End
###################################
gen_graphs = []

m=2
g=generate_g_2(m)
g_old = g.copy()
gen_graphs.append(g_old)
print gen_graphs

print '*******************************'

print "Is connected?", is_connected(g)

print "Is irreducible?", Irreducible(g)

reset_g(g,m)

print "Is hubbard_type?", Hubbard_diagram(g,m)


label_b_line(g,2)
print 'Labels = ', g.es["Label"] 

ami_input = AMI_Input(g,m)
print 'AMI_INPUT --> ', ami_input

print 'Is AMI result is zero?', AMI_zero(ami_input)
print '*******************************'
print


#print 'nodes tags =', g.vs["name1"]
#print 'visited = ', g.vs["visited"] 
#print 'spin = ', g.vs["Spin"]  
#print 'F_or_B = ',  g.es["F_or_B"]
#print 'INT or EXT = ', g.es["INT_or_EXT"] 
visual_style = {} 
color_dict = {"m": "black", "f": "white"}
visual_style["vertex_size"] = 20
visual_style["vertex_label"] = g.vs["name1"]
visual_style["edge_label"] = g.es["F_or_B"]
plot(g, '1.eps', **visual_style)  # Plot g


#start_time=time.time()
#b= AMI_arrays_out(ami_input,m)   # Generate S, P, and R arrays and store in b
#print("--- Construction_time ---", (time.time() - start_time))  # Compute the construction time
#print 'S = ', b[0]
#print 'P_freq = ', b[1]
#print 'P_mnta = ', b[2]
#print 'R_freq = ', b[3]
#print 'R_mnta = ', b[4]



remove_int_line(g)  # Remove one of the interaction lines
g_old = g.copy()
m=m-1

#print 'nodes tags =', g.vs["name1"]
#print 'visited = ', g.vs["visited"] 
#print 'spin = ', g.vs["Spin"]  
#print 'F_or_B = ',  g.es["F_or_B"]
#print 'INT or EXT = ', g.es["INT_or_EXT"] 

flag_iso=False
for i in range (0,len(gen_graphs)):
  if g_old.isomorphic_vf2(gen_graphs[i], return_mapping_12=False, return_mapping_21=False):
    graph_count=i
    flag_iso=True
    break
if flag_iso==False:
  gen_graphs.append(g_old)

#print flag_iso
#print gen_graphs  
    
print "Is connected?", is_connected(g)

print "Is irreducible?", Irreducible(g)

reset_g(g,m)

print "Is hubbard_type?", Hubbard_diagram(g,m)


print g
label_b_line(g,m)
#print 'Labels = ', g.es["Label"] 
ami_input = AMI_Input(g,m)


print 'AMI_INPUT --> ', ami_input

print 'Is AMI result is zero?', AMI_zero(ami_input)
print '*******************************'
print


#print 'visited = ', g.vs["visited"] 
#print 'spin = ', g.vs["Spin"]  
#print 'F_or_B = ',  g.es["F_or_B"]
#print 'INT or EXT = ', g.es["INT_or_EXT"] 


visual_style = {} 
color_dict = {"m": "black", "f": "white"}
visual_style["vertex_size"] = 20
visual_style["vertex_label"] = g.vs["name1"]
visual_style["edge_label"] = g.es["F_or_B"]
plot(g, '2.eps', **visual_style)  # Plot g

#start_time=time.time()
#b= AMI_arrays_out(ami_input,m)   # Generate S, P, and R arrays and store in b
#print("--- Construction_time ---", (time.time() - start_time))  # Compute the construction time
#print 'S = ', b[0]
#print 'P_freq = ', b[1]
#print 'P_mnta = ', b[2]
#print 'R_freq = ', b[3]
#print 'R_mnta = ', b[4]


add_int_line(g)  # Add an the interaction line
g_old = g.copy()
m=m+1

#print '****'
#print 'visited = ', g.vs["visited"] 
#print 'spin = ', g.vs["Spin"]  
#print 'F_or_B = ',  g.es["F_or_B"]
#print 'INT or EXT = ', g.es["INT_or_EXT"] 
#print "Is hubbardtype?", Hubbard_diagram(g,m)
flag_iso=False
for i in range (0,len(gen_graphs)):
  if g_old.isomorphic_vf2(gen_graphs[i], return_mapping_12=False, return_mapping_21=False):
    graph_count=i
    flag_iso=True
    break
if not flag_iso:
  gen_graphs.append(g)

#print flag_iso

print "Is connected?", is_connected(g)

print "Is irreducible?", Irreducible(g)

reset_g(g,m)

print "Is hubbard_type?", Hubbard_diagram(g,m)


#print g
label_b_line(g,m)
#reset_g(g,2)
ami_input = AMI_Input(g,m)
print 'AMI_INPUT --> ', ami_input
#print 'Labels = ', g.es["Label"] 

#print 'visited = ', g.vs["visited"] 
#print 'spin = ', g.vs["Spin"]  
#print 'F_or_B = ',  g.es["F_or_B"]
#print 'INT or EXT = ', g.es["INT_or_EXT"] 


visual_style = {} 
color_dict = {"m": "black", "f": "white"}
visual_style["vertex_size"] = 20
visual_style["vertex_label"] = g.vs["name1"]
visual_style["edge_label"] = g.es["F_or_B"]
plot(g, '3.eps', **visual_style)  # Plot g

#print gen_graphs

ami_input = AMI_Input(g,m)
print 'AMI_INPUT --> ', ami_input
print 'Is AMI result is zero?', AMI_zero(ami_input)
print '*******************************'
print

#start_time=time.time()
#b= AMI_arrays_out(ami_input,m)   # Generate S, P, and R arrays and store in b
#print("--- Construction_time ---", (time.time() - start_time))  # Compute the construction time
#print 'S = ', b[0]
#print 'P_freq = ', b[1]
#print 'P_mnta = ', b[2]
#print 'R_freq = ', b[3]
#print 'R_mnta = ', b[4]






