#!/usr/bin/env python
# coding: utf-8
## params ##

import numpy as np
import scipy.io as sio

global vertex_name_list, time_scale, c, r_default, epsilon, d_valid

simu_step = 1000
time_scale = 100
h = 1e3

c = 1
r_default = 0.8
r_new = 0.8
l_new = 1
epsilon = 0.1

###################################################################

handle = 'Q'            # Q - information entropy; M - missed probability; U - Li & Dahleh (2023)

solver = 1              # 1 - brute force; 2 - simulated annealing; 3 - sequential optimal

k_max = 3               # for iteration: maximum number of sensors to be placed

###################################################################

# Model: O2004 / OS2001 #

Matlab_file = sio.loadmat('data/O2004.mat')   # 'data/OS2001.mat'

vertex_name_list = list(Matlab_file['A_name'])
A = np.matrix(Matlab_file['A_full'])
d = np.array([0]*len(vertex_name_list))             # Data availability vector 

d = np.array([1,0,1,0,0,0,0])     # Existing data availability -- vertex_name_list = ['Br','B','P','Dr','D','Pd','A']

r = d*(1+r_default)-1                               # Data reliability vector 
l = np.array([100]*len(vertex_name_list))           # Data length vector 
d_valid = np.array(list(Matlab_file['A_v2_full'][0]))# 1: can be a sensor; 0: cannot be a sensor

l = l/time_scale
r = np.array(r,dtype=float)
l = np.array(l,dtype=float)

