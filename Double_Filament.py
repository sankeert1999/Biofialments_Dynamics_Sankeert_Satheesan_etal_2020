###################################################################################################################################################################
######## This program gives the output for the two-filament system trajectories.######################################################################################
# ########## Developed by Sankeert Satheesan (October 2020)####################################################################################################### 


import random
import random
import math
import numpy as np
import csv
w = 0
x = np.zeros((10000000,1))                           #array for saving the length time evolution for filament #1
y = np.zeros((10000000,1))                                                                                          #NUmber of datapoints is T_Stop - T_Trans
t_exp = np.zeros((10000000,1))                       #array for saving the iteration time
t_g = np.zeros((10000000,1))                         #array for saving the gillepsie time
#######################################################################################################################################################################
C = 0.0
L_1 = []                                        #array representing the state of the growing tubule #1
L_2 = []                                        #array representing the state of growing tubule #2
N = 1000                                        #Number of subunits (Free subunit pool)
l_1 = 0                                         #Initail filament #1
l_2 = 0                                         #Initial filament #2
NG = 0
t = 0.0                                         #Initializing the reaction time
dt = 0.01
tau = 0.0           
texp = 0.0                                      #Initializing the experimental time
T_trans = 0
tstop = 10000000                                #run time for the programme   
q_1 = 0.075                                     #Second order rate constant (for growth) if the monomer is a T
q_2 = 0.00075                                   #Second order rate constant (for growth) if the monomer is a D
G_T = 24.0                                      #Decay constant for GTP at end
G_D = 290.0                                     #Decay constant for GTP at end
h = 0.001                                         # rate of hydrolysis
H = 0.9
i = 0
m_1 = 0                                         # for randomly choosing the position of hydrolysis
m_2 = 0                                         # for randomly choosing the position of hydrolysis
a_1 = 0                                         #Initialization of rate constant
a_2 = 0
a_3 = 0
a_4 = 0
a_5 = 0
a_6 = 0
a_7 = 0
a_0 = 0      
a_sec = 0
a_thi = 0
a_fou = 0
a_fiv = 0
a_six = 0
for i in range(l_1):
    L_1.append(0)                               #initialization
for i in range(l_2):
    L_2.append(0)                               #initialization
#########################################################################################################################################################################
for texp in range(tstop):
    
##        Nf = float(N-NG)                      #having the cont of the number of free GTP monomers in the pool
    l1p = []                                    #array tracking the indexs with T monomers for filament_1
    l2p = []                                    #array tracking the indexs with T monomers for filament_2
    N_1 = L_1.count(0)                          #counting the number of GTP states (i.e. number of zeros)
    N_2 = L_2.count(0)
##        nfs = float(Nf-l)
#######################################################################################################################################################################
    
    a_2 = float(h*(N_1))                        # rate of hydrolysis       of Filament #1
    
    a_5 = float(h*(N_2))                        # rate of hydrolysis       of Filament #2 
    
    a_7 = float(H*(NG))                         #rate of neuclotide exchange

########################################################################################################################################################################    
    if l_1!=0:                                  #a nested if loop
        if L_1[-1] == 0:                        #if condition to select the decay constant according to the GTP/GDP states for filament#1
            a_3=G_T
            a_1 = float(q_1*(N - NG -(l_1 + l_2)))        #rate updation with the desired fromulae of particular length regulation  of Filament #1 if monomer is T
        elif L_1[-1] == 1:
            a_3=G_D
            a_1 = float(q_2*(N - NG -(l_1 + l_2)))        #rate updation with the desired fromulae of particular length regulation  of Filament #1 if monomer is D
    else:
        a_2 = 0
        a_3 = 0
        a_1 = float(q_1*(N - NG -(l_1 + l_2)))        #rate updation with the desired fromulae of particular length regulation  of Filament #1 if monomer is T
##########################################################################################################################################################################
    if l_2 != 0:
        if L_2[-1] == 0:                        #if condition to select the decay constant according to the GTP/GDP states for filament #2 
            a_6=G_T
            a_4 = float(q_1*(N - NG -(l_1 + l_2)))        #rate updation with the desired fromulae of particular length regulation  of Filament #2 if monmer is T
        elif L_2[-1] == 1:
            a_6=G_D
            a_4 = float(q_2*(N - NG -(l_1 + l_2)))        #rate updation with the desired fromulae of particular length regulation  of Filament #2 if monomer is D
    else :
        a_5 = 0
        a_6 = 0
        a_4 = float(q_1*(N - NG -(l_1 + l_2)))        #rate updation with the desired fromulae of particular length regulation  of Filament #2 if monomer is T
#################################################################################################################################################################################
    a_sec = a_1 + a_2
    a_thi = a_1 + a_2 + a_3
    a_fou = a_1 + a_2 + a_3 + a_4
    a_fiv = a_1 + a_2 + a_3 + a_4 + a_5
    a_six = a_1 + a_2 + a_3 + a_4 + a_5 + a_6
    a_0 = a_1 + a_2 + a_3 + a_4 + a_5 + a_6 + a_7
    tam = random.uniform(0,1)
    tau = (1.0/a_0)*math.log(1.0/tam)           #Generating a random reaction time gillepsie formalism
    t = t + tau                                 # updation of the reaction time
    k = random.uniform(0,1)                     #picking the random number for choosing the reaction
##    print(R,G)
################################################################################################################################
    for j in range(0,(len(L_1))):
        if L_1[j] == 0:                         # appending all the indexs with T monomers in tracker array for filament#1
            l1p.append(j)
    for u in range(0,(len(L_2))):
        if L_2[u] == 0:                         # appending all the indexs with T monomers in tracker array for filament#1
            l2p.append(u) 
##########################################################################################################################    
    if 0<=k*a_0 and k*a_0<=a_1:                 #for the reaction 1 (hydrolysis)
        l_1 = l_1 + 1                           
        L_1.insert(l_1,0)
        
    elif a_1<=k*a_0 and k*a_0<=a_sec:           #Condition for the occurence of growth reaction 
        m = random.randint(0,(len(l1p) - 1))    #picking the random position
        p1 = l1p[m]                             #saving in a dummy variable         For Filament#1
        L_1[p1] = 1                             #updating the main lattice

    elif a_sec<=k*a_0 and k*a_0<=a_thi:         #Condition for the occurence of decay reaction
        l_1 = l_1 - 1
        L_1.pop(-1)
        if a_3 == G_D:                          #Introduction of the if loop as the number of GDP should be updated only in condition with the end being G_D
            NG = NG + 1                         #Updating the number of GDP Subunits
########################################################################################################################            
    elif a_thi<=k*a_0 and k*a_0<=a_fou:         #for the reaction 1 (hydrolysis)
        l_2 = l_2 + 1                           
        L_2.insert(l_2,0)
        
    elif a_fou<=k*a_0 and k*a_0<=a_fiv:         #Condition for the occurence of growth reaction 
        m = random.randint(0,(len(l2p) - 1))    #picking the random position
        p2 = l2p[m]                             #saving in a dummy variable         For filament#2
        L_2[p2] = 1                             #updating the main lattice
    elif a_fiv<=k*a_0 and k*a_0<=a_six:         #Condition for the occurence of decay reaction
        l_2 = l_2 - 1
        L_2.pop(-1)
        if a_6 == G_D:                          #Introduction of the if loop as the number of GDP should be updated only in condition with the end being G_D
            NG = NG + 1                         #Updating the number of GDP Subunits            
    elif a_six<=k*a_0 and k*a_0<=a_0:           #conditon for the occurence of a dehyrolysis reaction
        NG = NG - 1

    
###############################################################################################################################3
#    if  l!=0:                                  #Only noting down values for L not equal to zero 
        
##
    if texp>=T_trans and texp<=tstop-dt:
        x[w] = l_1
        y[w] = l_2
        t_exp[w] = texp
        t_g[w] = t            
        w = w+1
  
#######################################################################################################################################################################

np.savetxt("h=0.001.csv", np.column_stack((t_exp,t_g,x,y))) 




    

    
    
