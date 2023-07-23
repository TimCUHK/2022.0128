#!/usr/bin/env python
# coding: utf-8
## Functions ##

import networkx as nx
import operator
import math
import numpy as np
import random
from itertools import permutations 
from itertools import combinations 
import copy
import scipy.io as sio
import time

from params import *

vertex_name_list = list(Matlab_file['A_name'])
A = np.matrix(Matlab_file['A_full'])
d = np.array([0]*len(vertex_name_list))             # Data availability vector 
r = d*(1+r_default)-1                               # Data reliability vector 
l = np.array([100]*len(vertex_name_list))           # Data length vector 
d_valid = np.array(list(Matlab_file['A_v2_full'][0]))# 1: can be a sensor; 0: cannot be a sensor


l = l/time_scale
r = np.array(r,dtype=float)
l = np.array(l,dtype=float)

# Function: Initialization # 
def Initialization(A,d,r,l,info_read):
    
    G = nx.Graph()

    for i in range(len(d)):
        G.add_node(i, Name = vertex_name_list[i]) 

    AT = A.transpose()

    if len(AT) != len(d):
        print('Error: Matrix and data availability vector length do not match.')

    # Parameter Identification    
    para = []
    for i in range(len(d)):

        if np.sum(AT[i,:]) == 0:      # Node has no predecessor
            G.add_node(i, Parameter = 'True')
            para.append(i) 
        else:
            G.add_node(i, Parameter = 'False')

    # Reliability Vector Calculation

    r,FVS = Calculation_R(A,d,r)
    if info_read == 1:
        print('[Data Reliability Vector - r]')
        print([round(x, 1) for x in r])

        print('[FVS]')
        print(FVS)

    return r,FVS,para,G

# Function: Calculate Data_Availability_Partition (Oliva 2004) #

def Data_Availability_Partition(AT,d,r,node_ind,G):
    
    # Adjacency Matrix, Data availability vector, Data reliability vector, Node
    y = node_ind
    x = []
    beta = []
    eq = []
    
    # initial coloring -- necessary; update for each use of function
    for i in range(len(G)):
        G.nodes[i]['color'] = 'white'
        if i == node_ind:
            G.nodes[i]['color'] = 'black'
    
    Q_list = []
    Q_list.append(node_ind)
   
    while len(Q_list)>0:
        
        head = Q_list[0]
        # Find the predecessor of head
        seq = [k for k, e in enumerate(AT.tolist()[head]) if e != 0]
        
        for i in seq:
            
            if G.nodes[i]['color'] == 'white':
                
                G.nodes[i]['color'] = 'black'  # explored
                
                if d[i] == 1:   # data available
                    x.append(i)
                elif G.nodes[i]['Parameter'] == 'True':
                    beta.append(i)
                else:
                    eq.append(i)
                    Q_list.append(i)
    
        Q_list = Q_list[1:]     
  
        G.nodes[head]['color'] = 'black'  # explored
    
    return y,x,beta,eq
            # dependent variable, data-available nodes, parameters, intermediate nodes

# Objective Function # -- information entropy (PL2012)

def Objective_Func_Q(A,d,r,l,FVS,para,G):
    
    # L #
    d_indices = [i for i, x in enumerate(d) if x == 1]
    L = np.zeros((len(d_indices),len(d)))
    for i in range(len(d_indices)):
        L.itemset((i, d_indices[i]), 1)

    # phi #
    phi_index = d_indices + para + FVS
    phi_index = list(dict.fromkeys(phi_index))
    phi = np.zeros((len(d),len(phi_index)))

    AT = A.transpose()
    for i in range(len(d)):

        x_i,x_d_i,alpha_i,beta_i = Data_Availability_Partition(AT,d,r,i,G)  
             # dependent variable, data-available nodes, parameters, intermediate nodes
        todo = x_d_i + alpha_i + list(set(beta_i) & set(FVS)) # data nodes, parameters, FVS
        todo = list(dict.fromkeys(todo))
        for j in range(len(todo)):
            phi.itemset((i, phi_index.index(todo[j])), 1)  

    # Sigma #
    Sigma = np.diag(np.float_power(r, -1))

    # Q # 
    Lphi = np.dot(L,phi)
    core = np.linalg.inv(np.linalg.multi_dot([L,Sigma,L.transpose()]))             
    Q = np.linalg.multi_dot([Lphi.transpose(),core,Lphi])

    # Eigenvalues #

    w,v = np.linalg.eig(Q)
    non_indices = [i for i, x in enumerate(w) if abs(x) > 1e-6]

    value = np.real(np.prod([w[i] for i in non_indices]))

    return Q,value

# Objective Function # -- miss probability (DC2003)

def Objective_Func_M(A,d,r,l,FVS,para,G):
    
    Mv = []
    AT = A.transpose()
    
    for i in range(len(d)):

        x_i,x_d_i,alpha_i,beta_i = Data_Availability_Partition(AT,d,r,i,G)  

        todo = x_d_i + alpha_i + list(set(beta_i) & set(FVS)) # data nodes, parameters, FVS
        todo = list(dict.fromkeys(todo))
        if len(todo) == 0:
            m = 0
        else:    
            m = (1-r[i])*len(alpha_i)/(len(todo))

        Mv.append(m)

    value = -sum(Mv)

    return Mv,value

# Objective Function # -- individual utility (LD2023)

def Objective_Func_U(A,d,r,l,FVS,para,G):
    
    AT = A.transpose()
    len_vertex = len(d)
    
    value = 0
    best_seq = []
    
    # Permutation of all sequence --> Only one data-available variable is removed from d at each time.
    sequence = list(permutations([ee for ee, e in enumerate(d) if e != 0]))

    for i in range(len(sequence)):
        
        temp_sequence = sequence[i]
        
        # Copy initial input
        dd = list(np.copy(d))
        rr = list(np.copy(r))
        ll = list(np.copy(l))
        AAT = np.copy(AT)
        
        temp_value = 0
        
        reference_vector = list(range(0, len_vertex))
        
        for j in range(len(temp_sequence)):
            
            floating_index = reference_vector.index(temp_sequence[j])
            
            temp_utility,temp_DAP_info = Calculation_Utility(AAT,dd,rr,ll,floating_index,G)
            
            temp_value = temp_value + temp_utility
            
            y = temp_DAP_info[0]
            x = temp_DAP_info[1]
            beta = temp_DAP_info[2]
            eq = temp_DAP_info[3]
            
            # Reshaping Matrix and Vectors
            del_col = [y] + beta + eq   # Entries in A to be deleted
                        
            AAT = np.delete(AAT,del_col,1)
            AAT = np.delete(AAT,del_col,0)
            
            dd = np.delete(dd,del_col,0)
            rr = np.delete(rr,del_col,0)
            ll = np.delete(ll,del_col,0)

            reference_vector = list(np.delete(reference_vector,del_col))
        
        # Update U
        if temp_value > value:
            value = temp_value
            best_seq = temp_sequence
    
    return best_seq,value

# Function: Calculate Reliability Score #

def Calculation_R(A,d0,r0):
    
    global r_default
    
    d = d0.copy()
    r = r0.copy()
    
    AT = A.transpose()
    GG = nx.from_numpy_matrix(A,create_using=nx.DiGraph)
    
    # -- FVS Identification
    FVS = FVS_ID(GG,d)

    if len(FVS)!= 0:
        for i in range(len(FVS)):
            if r[FVS[i]] == -1:
                r[FVS[i]] = r_default
                d[FVS[i]] = 1   ##########!!!!!!! 
    
    # -- Assignment
    while True:
       
        if sum([k == -1 for k in r]) == 0: 
            break

        temp_node = random.choice([k for k in range(len(d))]) #(list(G.nodes()))

        if r[temp_node] != -1: 
            continue

        predecessor = [ee for ee, e in enumerate(AT[temp_node,:].tolist()[0]) if e != 0]

        
        continue_flag = 0
        if len(predecessor) == 0:
            new_r = r_default
        else:
            if -1 in [r[k] for k in predecessor]:
                continue_flag = 1 
            new_r = np.prod([r[k] for k in predecessor])
        
        if continue_flag == 1:
            continue
            
        if r[temp_node] != new_r:
            r[temp_node] = new_r
            # reset
            successor = [ee for ee, e in enumerate(A[temp_node,:].tolist()[0]) if e != 0]
            if len(successor) != 0:
                for k in successor:
                    if d[k]!= 1:
                        r[k] = -1
    
    return r,FVS

# Function: Identifying Feedback Vertex Set #

def FVS_ID(GG,d):
   
    FVS = []
    cycle = list(nx.simple_cycles(GG))
    cycle.sort(key=len,reverse = True)
    for i in range(len(cycle)):
        
        if sum(d[k] for k in cycle[i]) != 0:
            continue
            
        if len(set(FVS) & set(cycle[i])) == 0:
            max_d = 0
            max_d_ind = []
            for j in range(len(cycle[i])):
                if GG.degree(cycle[i][j]) > max_d:
                    max_d = GG.degree(cycle[i][j])
                    max_d_ind = j
                    
            FVS.append(cycle[i][max_d_ind])
            
    return FVS

# Function: Calculate Transitive Closure of a matrix (2*2 list) #

def Warshall(a):
    assert (len(row) == len(a) for row in a)
    n = len(a)
    for k in range(n):
        for i in range(n):
            for j in range(n):
                a[i][j] = a[i][j] or (a[i][k] and a[k][j])
    return a

# Function: Calculate New Data Utility (NDU) #

def Calculation_NDU(index,r_new,l_new,AT,d,r,l,G):
    
    dd = list(np.copy(d))
    dd[index] = 1
    
    rr = list(np.copy(r))
    rr[index] = r_new
    
    ll = list(np.copy(l))
    ll[index] = l_new
    
    rr,FVS = Calculation_R(AT,dd,rr)
        
    NDU = Calculation_DAVV(AT,dd,rr,ll,G) - Calculation_DAVV(AT,d,r,l,G)

    return NDU 
     
# Function: Calculate Data Ability Vector Value (DAVV) #

def Calculation_DAVV(AT,d,r,l,G):
    
    len_vertex = len(d)
    
    DAVV = 0
    best_seq = []
    
    # Permutation of all sequence --> Because only one available data is removed from d at each time.
    sequence = list(permutations([ee for ee, e in enumerate(d) if e != 0]))

    for i in range(len(sequence)):
        
        temp_sequence = sequence[i]
        
        # Copy initial input
        dd = list(np.copy(d))
        rr = list(np.copy(r))
        ll = list(np.copy(l))
        AAT = np.copy(AT)
        
        temp_DAVV = 0
        
        reference_vector = list(range(0, len_vertex))
        
        for j in range(len(temp_sequence)):
            
            floating_index = reference_vector.index(temp_sequence[j])
            
            temp_utility,temp_DAP_info = Calculation_Utility(AAT,dd,rr,ll,floating_index,G)
            
            temp_DAVV = temp_DAVV + temp_utility
            
            y,x,beta,eq = Data_Availability_Partition(AAT,dd,rr,floating_index,G)
            
            # Reshaping Matrix and Vectors
            del_col = [y] + beta + eq   # Entries in A to be deleted
                        
            AAT = np.delete(AAT,del_col,1)
            AAT = np.delete(AAT,del_col,0)
            
            dd = np.delete(dd,del_col,0)
            rr = np.delete(rr,del_col,0)
            ll = np.delete(ll,del_col,0)


            reference_vector = list(np.delete(reference_vector,del_col))
            
        # Update DAVV
        if temp_DAVV > DAVV:
            DAVV = temp_DAVV
            best_seq = temp_sequence

    return DAVV

# Function: Select the Best N Variables for New Data Acquisition (OptN) #

def Calculation_OptN(AT,d,r,l,num_ND,r_new,l_new,G):
    
    len_vertex = len(d)
    
    # Combination of all slots --> Because only one available data is added to d at each time.
    if len([ee for ee, e in enumerate(d) if e == 0]) < num_ND:
        print('Error! Less than N variables have no data.')
        return [0]
    
    sequence = list(combinations([ee for ee, e in enumerate(d) if e == 0],num_ND))
    
    OptN = []
    OptN_U = 0
    
    for i in range(len(sequence)):

        temp_sequence = sequence[i]
        
        # Copy initial input
        dd = list(np.copy(d))
        rr = list(np.copy(r))
        ll = list(np.copy(l))
        AAT = np.copy(AT)
        
        temp_OptN_U = 0
               
        for j in range(len(temp_sequence)):
             
            floating_index = temp_sequence[j]   # no reduction of matrix, no need of reference vector for shrinkage
            
            temp_NDU = Calculation_NDU(floating_index,r_new,l_new,AAT,dd,rr,ll,G)
        
            temp_OptN_U = temp_OptN_U + temp_NDU
            
            AAT = AAT
            dd[floating_index] = 1
            rr[floating_index] = r_new
            ll[floating_index] = l_new
        
        # Update OptN
        if temp_OptN_U == OptN_U:
            OptN.append(temp_sequence)  # more than one possible sequence 
            
        if temp_OptN_U > OptN_U:
            OptN_U = temp_OptN_U
            OptN = list([temp_sequence])

    return OptN   # index in original d vector
            
# Function: Select the Optimal Sequence for New Data Acquisition (OptSeq) #

def Calculation_OptSeq(AT,d,r,l,r_new,l_new,G):
    
    len_vertex = len(d)
    
    remain_seq = [ee for ee, e in enumerate(d) if e == 0]
    
    OptSeq = []
    
    # Copy initial input  
    dd = list(np.copy(d))
    rr = list(np.copy(r))
    ll = list(np.copy(l))
    AAT = np.copy(AT)
        
    for i in range(len(remain_seq)):
         
        temp_OptN = Calculation_OptN(AAT,dd,rr,ll,1,r_new,l_new,G)
        if temp_OptN == []:
            break
        
        ind = temp_OptN[0][0]
        
        AAT = AAT
        dd[ind] = 1
        rr[ind] = r_new
        ll[ind] = l_new
        
        OptSeq.append(ind)
    
    return OptSeq

# Function: Calculate Utility for node that has available data #

def Calculation_Utility(AT,d,r,l,node_ind,G):
    
    global c,epsilon
    
    y,x,beta,eq = Data_Availability_Partition(AT,d,r,node_ind,G)
    
    if len(beta) == 0:    
        if len(eq) != 0:    # Not a parameter. Cannot do calibration.
            temp = 0
        else:               # Is a parameter
            temp = r[y] * l[y]*epsilon
    else:
       
        temp = (len(eq) + 1)**c/len(beta) * min([r[k]*l[k] for k in x+[y]])
        
        temp = (len(eq) + 1)**c/len(beta) * min([l[k] for k in x+[y]])*sum([r[k] for k in x+eq+[y]])/(len(eq)+1)
        
    Utility = temp
    DAP_info = [y,x,beta,eq]
    
    return Utility, DAP_info

# Function: Output Objective values, Given a candidate of d-r-l #

def Func_Output_single_d(A,d0,r0,l0):
    global c,r_default
    
    d = d0.copy()
    r = r0.copy()
    l = l0.copy()
    # Initialization #
    print('[Data Availability Vector - d]')
    print(d)
    print('[Data Length Vector - l]')
    print([round(x,1) for x in l])
    r,FVS,para,G = Initialization(A,d,r,l,1)
    # Objective Function #
    Q, value1 = Objective_Func_Q(A,d,r,l,FVS,para,G)
    M, value2 = Objective_Func_M(A,d,r,l,FVS,para,G)
    best_seq, value3 = Objective_Func_U(A,d,r,l,FVS,para,G)
    print('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*')
    print('[Objective Value - Q]')
    print(round(value1, 2))
    print('[Objective Value - M]')
    print(round(value2, 2))
    print('[Objective Value - U]')
    print(round(value3, 7))
    
    return value1,value2,value3

# Solver: Simulated Annealing search the optimal d, given k and d0 #

def Func_Simulated_Annealing_Optimal(A,d0,r0,l0,k,r_new,l_new,handle,simu_step,info_read):

    global c,r_default,h,d_valid
    
    if handle == 'Q':
        obj_f = Objective_Func_Q
    if handle == 'M':
        obj_f = Objective_Func_M
    if handle == 'U':
        obj_f = Objective_Func_U
        
    candidate_indices = [i for i, x in enumerate(d0) if (x == 0) & d_valid[i] == 1]
    if len(candidate_indices) < k:
        print('Warning! k is too large.')
        return 0
    
    pool = list(combinations(candidate_indices, k))
    visited = [0]*len(pool)

    # initial value #
    d = d0.copy()
    r = r0.copy()
    l = l0.copy()
    r = np.array(r,dtype=float)
    l = np.array(l,dtype=float)
    
    ini = random.randint(0, len(pool)-1)
    
    for j in range(len(pool[ini])):
        d[pool[ini][j]] = 1
        r[pool[ini][j]] = r_new
        l[pool[ini][j]] = l_new

    # Initialization #
    r,FVS,para,G = Initialization(A,d,r,l,0)
    # Objective Function #
    out1, value = obj_f(A,d,r,l,FVS,para,G)
    
    minimum_obj = -value
    current = list(pool[ini])   # current placement
    d_optimal = d
    
    t_max = simu_step
    t = 0
    
    while t < t_max:
        t = t+1
        
        trial_time = 0
        while trial_time < 100:
            trial_time = trial_time + 1
            trial = random.randint(0, len(d)-1)
            if trial in candidate_indices and trial not in current:
                break
                
        current[random.randint(0, len(current)-1)] = trial
        
        d = d0.copy()
        r = r0.copy()
        l = l0.copy()
        r = np.array(r,dtype=float)
        l = np.array(l,dtype=float)
        
        for j in range(len(current)):
            d[current[j]] = 1
            r[current[j]] = r_new
     
        # Initialization #
        r,FVS,para,G = Initialization(A,d,r,l,0)
        # Objective Function #
        out1, value = obj_f(A,d,r,l,FVS,para,G)
    
        # Judging condition #
        delta_value = math.exp((minimum_obj + value)/(t*h))
        
        if delta_value > random.uniform(0, 1):
            minimum_obj = -value
            d_optimal = d.copy()
            if (info_read == 1) & (t%100 == 0):    
                print('[FSAO - Iteration '+ str(t)+'/'+str(t_max)+']:' + str(d) + '--' + str(minimum_obj))
            
    return d_optimal

# Solver: Sequentially optimal d, given k and d0 #

def Func_Sequential_Optimal(A,d0,r0,l0,k,r_new,l_new,handle,info_read):
    global c,r_default,d_valid
    
    candidate_indices = [i for i, x in enumerate(d0) if (x == 0) & d_valid[i] == 1]
    if len(candidate_indices) < k:
        print('Warning! k is too large.')
        return 0
    
    d = d0.copy()
    r = r0.copy()
    l = l0.copy()
    r = np.array(r,dtype=float)
    l = np.array(l,dtype=float)
    
    for i in range(k):
        
        round_result = Func_Brute_Force_Optimal(A,d,r,l,1,r_new,l_new,handle,info_read)
        
        d[round_result[0]] = 1
        r[round_result[0]] = r_new
        l[round_result[0]] = l_new
    
        if info_read == 1:
            print('[FSO - Iteration '+ str(i+1)+'/'+str(k)+']:' + str(round_result[0]))
        
    return d

# Solver: Brute force search the optimal d, given k and d0 #

def Func_Brute_Force_Optimal(A,d0,r0,l0,k,r_new,l_new,handle,info_read):
    global c,r_default,d_valid
    
    if handle == 'Q':
        obj_f = Objective_Func_Q
    if handle == 'M':
        obj_f = Objective_Func_M
    if handle == 'U':
        obj_f = Objective_Func_U
        
    candidate_indices = [i for i, x in enumerate(d0) if (x == 0) & d_valid[i] == 1]
    if len(candidate_indices) < k:
        print('Warning! k is too large.')
        return 0
    
    pool = list(combinations(candidate_indices, k))
    
    maximum_obj = -12345
    maximum_d = []
    
    for i in range(len(pool)):
        
        d = d0.copy()
        r = r0.copy()
        l = l0.copy()
        r = np.array(r,dtype=float)
        l = np.array(l,dtype=float)
        for j in range(len(pool[i])):
            d[pool[i][j]] = 1
            r[pool[i][j]] = r_new
            l[pool[i][j]] = l_new

        # Initialization #
        r,FVS,para,G = Initialization(A,d,r,l,0)

        # Objective Function #
        out1, value = obj_f(A,d,r,l,FVS,para,G)
        
        if abs(value - maximum_obj)< 1e-3:
            maximum_d.append([x for x in pool[i]])
        elif value > maximum_obj:
            maximum_obj = value
            maximum_d = [[x for x in pool[i]]]
        if (info_read == 1) & (i%100 == 0) & (i>99):    
            print('[FBFO - Iteration '+ str(i+1)+'/'+str(len(pool))+']:' + str([x for x in pool[i]]) + '--' + str(value))

    return maximum_d 


