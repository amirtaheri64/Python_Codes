from Symbolic_multi_AMI_new import arr_comp_new

def AMI_INPUT_DEC(AMI_INPUT):
  a = [] 
  num = []
  result_freq = []
  #ovl_sgn = 1
  for i in range(0,len(AMI_INPUT)):
    a.append(0)
    num.append(0)
    result_freq.append(0)
    
  l = 0 
  for i in range (0,len(AMI_INPUT)):
    k = 0
    
    if a[i] == 0:
      a[i] = 1
      
      result_freq[l] = AMI_INPUT[i]  
      for j in range (i+1, len(AMI_INPUT)):
        comp_freq = arr_comp_new(AMI_INPUT[i],AMI_INPUT[j])
        if comp_freq[0]:
          #ovl_sgn = ovl_sgn*comp_freq[1]
          a[j] = 1
          k = k+1
      num[l] = k+1
      #print result[l], num[l], l
      l = l + 1
 
  return result_freq, num#, ovl_sgn

def AMI_zero(ami_input):
  ami_input_dec = AMI_INPUT_DEC(ami_input)
  #print ami_input_dec

  l = 0
  for i in range (0, len(ami_input_dec[1])):
    if ami_input_dec[1][i]!=0:
      #print i, dec[1][i], dec[0][i]
      l = l + 1
    #print l
    dec_new = [[None]*l for i in range(0,2)]
    #print 'dec_new = ', dec_new
  for i in range(0,l):
    dec_new[0][i] = ami_input_dec[0][i]
    dec_new[1][i] = ami_input_dec[1][i] 
  #overal_sign = ami_input_dec[2]

  #print dec_new

  arr_rep = dec_new[0]
  multiplicity=dec_new[1]
  #print 'arr_rep = ', arr_rep
  #print 'multiplicity = ', multiplicity

  a=[]
  #for i in range (0,len(multiplicity)):
    #if multiplicity[i]==1:
      #flag=1
  flag=1
  for i in range (0,len(arr_rep)):
    if multiplicity[i]>1:
      
      for elements in range(0,len(arr_rep[i])-1):
        flag_new=1
        if arr_rep[i][elements]!=0:
          flag_new=0
          for j in range(0,len(arr_rep)):
            if i!=j:
              if arr_rep[j][elements]!=0:
                flag_new=1
        flag=flag*flag_new
        if flag==0:
          break
      if flag==0:
        break
  
    
  return flag#,a

labels=[ [1,0,0],[1,0,0],[0,1,0]]
print AMI_zero(labels)


  



