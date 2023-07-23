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

handle = 'U'            # Q - information entropy; M - missed probability; U - Li & Dahleh (2023)

k = 1                   # number of sensors to be placed

solver = 1              # 1 - brute force; 2 - simulated annealing; 3 - sequential optimal


# Model: O2004 / OS2001 #

Matlab_file = sio.loadmat('../data/O2004.mat')   # 'data/OS2001.mat'


