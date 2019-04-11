from Label_Pol import *
import time

###########################################################
#########################Main##############################
###########################################################   

m = 5  # Order of the diagram
step=1
#n_vg = N_V(m)  # Number of the nodes
#lg = L_F(m) + L_I(m) + 2  # Total number of the edges
start_time = time.time()
#file = open("start_from_1", "w")
for i in range (0,step):
  #txt ='' 
  g = generate(m)  # Define a graph of order m
  result_s = label_sys(g,m) # Generate the labels
  node_seq = result_s[1]
  #print node_seq
  #for j in range(0, m+1):
    #txt = txt + str(node_seq[j]) + '\t'
  #txt = txt + '\n' 
  #file.write(txt)
  #print
end_time = time.time()
#file.close()
print "Sytematic Time = ", (end_time-start_time)/step, "sec"
max_v=400   # Maximum number of random nodes in each try
start_time = time.time()
for i in range (0,step):
  g = generate(m)  # Define a graph of order m
  result_r = label_ran(g,m,max_v) # Generate the labels
end_time = time.time()
print "Random Time = ", (end_time-start_time)/step, "sec"
start_time = time.time()
#file = open("start_from_1", "w")
for i in range (0,step):
  #txt ='' 
  g = generate(m)  # Define a graph of order m
  result_f = label_f_line(g,m) # Generate the labels
  node_seq = result_f[1]
  #print node_seq
  #for j in range(0, m+1):
    #txt = txt + str(node_seq[j]) + '\t'
  #txt = txt + '\n' 
  #file.write(txt)
  #print
end_time = time.time()
#file.close()
print "f-line Time = ", (end_time-start_time)/step, "sec"

##############################
######Systematic Result#######
##############################
#Begin
print
print "****************************"
print "******Systematic Result*****"
print "****************************"
print "Label = ", result_s[0].es["Label"]
print "F_or_B = ", result_s[0].es["name2"]

#End
##############################

##############################
########Random Result#########
##############################
#Begin
print
print "****************************"
print "*******Random Result********"
print "****************************"
print "try_label = ", result_r[1]
print "Label = ", result_r[0].es["Label"]
print "F_or_B = ", result_r[0].es["name2"]

#End
##############################


##############################
########f_line Result#########
##############################
#Begin
print
print "****************************"
print "********f-line Result*******"
print "****************************"
print "Label = ", result_f[0].es["Label"]
print "F_or_B = ", result_f[0].es["name2"]

#End
##############################

#############################
############Plot##############
##############################
#Begin
visual_style = {}
color_dict = {"m": "black", "f": "white"}
visual_style["vertex_size"] = 20
visual_style["vertex_label"] = result_s[0].vs["name1"]
visual_style["edge_label"] = result_s[0].es["name2"]
plot(g, **visual_style)
#End
###############################

#End Main
###########################################################


