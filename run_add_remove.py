##################################
##########Import Libraries########
##################################
#Begin
from add_remove import *
from is_connected import *
#End
###################################
  
m=2  # Choose the order to be zero
g = generate(m)   # Generate the zeroth order diagram

print '*********First***********'
visual_style = {} 
color_dict = {"m": "black", "f": "white"}
visual_style["vertex_size"] = 20
visual_style["vertex_label"] = g.vs["name1"]
visual_style["edge_label"] = g.es["F_or_B"]
plot(g,"0.eps", **visual_style)  # Plot g
plot(g, **visual_style)  # Plot g


add_int_line(g)
add_bubble(g)
visual_style["vertex_size"] = 20
visual_style["vertex_label"] = g.vs["name1"]
visual_style["edge_label"] = g.es["F_or_B"]
plot(g,"0.eps", **visual_style)  # Plot g
plot(g, **visual_style)  # Plot g
print 'g =', g
result_s=label_sys(g,3)
print result_s[0].es["Label"]
print bubble_finder(g)

print '*********Second*********'
visual_style = {} 
color_dict = {"m": "black", "f": "white"}
visual_style["vertex_size"] = 20
visual_style["vertex_label"] = g.vs["name1"]
visual_style["edge_label"] = g.es["F_or_B"]
plot(g,"1.eps", **visual_style)  # Plot g
plot(g, **visual_style)  # Plot g

print '************Third*********'
g_old = g.copy()
remove_int_line(g)  # Remove one of the interaction lines
visual_style = {} 
color_dict = {"m": "black", "f": "white"}
visual_style["vertex_size"] = 20
visual_style["vertex_label"] = g.vs["name1"]
visual_style["edge_label"] = g.es["F_or_B"]
plot(g,"2.eps", **visual_style)  # Plot g
plot(g, **visual_style)  # Plot g
conn = is_connected(g)
print conn


if conn==False:
  g = g_old.copy()
  print bubble_finder(g)

  visual_style = {} 
  color_dict = {"m": "black", "f": "white"}
  visual_style["vertex_size"] = 20
  visual_style["vertex_label"] = g.vs["name1"]
  visual_style["edge_label"] = g.es["F_or_B"]
  plot(g,"2_old.eps", **visual_style)  # Plot g
  plot(g, **visual_style)  # Plot g

print '**********Fourth*********'
g_old = g.copy()
remove_int_line(g)  # Remove one of the interaction lines
visual_style = {} 
color_dict = {"m": "black", "f": "white"}
visual_style["vertex_size"] = 20
visual_style["vertex_label"] = g.vs["name1"]
visual_style["edge_label"] = g.es["F_or_B"]
plot(g,"3.eps", **visual_style)  # Plot g
plot(g, **visual_style)  # Plot g
conn = is_connected(g)
print conn


if conn==False:
  g = g_old.copy()


  visual_style = {} 
  color_dict = {"m": "black", "f": "white"}
  visual_style["vertex_size"] = 20
  visual_style["vertex_label"] = g.vs["name1"]
  visual_style["edge_label"] = g.es["F_or_B"]
  plot(g,"3_old.eps", **visual_style)  # Plot g
  plot(g, **visual_style)  # Plot g
  

add_bubble(g)

add_int_line(g)

print '**********Fifth*********'
g_old = g.copy()
print bubble_finder(g)
visual_style = {} 
color_dict = {"m": "black", "f": "white"}
visual_style["vertex_size"] = 20
visual_style["vertex_label"] = g.vs["name1"]
visual_style["edge_label"] = g.es["F_or_B"]
plot(g,"4.eps", **visual_style)  # Plot g
plot(g, **visual_style)  # Plot g
conn = is_connected(g)
print conn





print '**********Sixth**********'
g_old = g.copy()
remove_int_line(g)  # Remove one of the interaction lines
visual_style = {} 
color_dict = {"m": "black", "f": "white"}
visual_style["vertex_size"] = 20
visual_style["vertex_label"] = g.vs["name1"]
visual_style["edge_label"] = g.es["F_or_B"]
plot(g, "5.eps", **visual_style)  # Plot g
plot(g, **visual_style)  # Plot g
conn = is_connected(g)
print conn

if conn==False:
  g = g_old.copy()


  #print g
  #print 'bubbles = ', bubble_finder(g)
  visual_style = {} 
  color_dict = {"m": "black", "f": "white"}
  visual_style["vertex_size"] = 20
  visual_style["vertex_label"] = g.vs["name1"]
  visual_style["edge_label"] = g.es["F_or_B"]
  plot(g, "5_old.eps",**visual_style)  # Plot g
  plot(g, **visual_style)  # Plot g

print '**********Seventh**********'
remove_bubble(g)  # Remove one of the bubbles
visual_style = {} 
color_dict = {"m": "black", "f": "white"}
visual_style["vertex_size"] = 20
visual_style["vertex_label"] = g.vs["name1"]
visual_style["edge_label"] = g.es["F_or_B"]
plot(g, "6.eps", **visual_style)  # Plot g
plot(g, **visual_style)  # Plot g
conn = is_connected(g)
print conn

if conn==False:
  g = g_old.copy()
  

  #print g
  #print 'bubbles = ', bubble_finder(g)
  visual_style = {} 
  color_dict = {"m": "black", "f": "white"}
  visual_style["vertex_size"] = 20
  visual_style["vertex_label"] = g.vs["name1"]
  visual_style["edge_label"] = g.es["F_or_B"]
  plot(g, "6_old.eps",**visual_style)  # Plot g
  plot(g, **visual_style)  # Plot g



