###############################################################################
#				Run Code for Mixed Membership Block Model
#				Authors: Alex Burnap, Efren Cruz, Xin Rong, Brian Segal
#				Date: April 1, 2013
###############################################################################

#---------------------------- Dependencies -----------------------------------#
import MMSB_model as model
import pymc
import numpy as np
import matplotlib.pyplot as plt

#---------------------------- Run Time Params --------------------------------#

# Probably going to try and use - flag conventions with __init__(self, *args, **kwargs)


#---------------------------- Load Data --------------------------------------#
data_matrix=np.loadtxt("../data/Y_alpha0.1_K5_N20.txt", delimiter=',')
num_people = 20
num_groups = 5
alpha = np.ones(num_groups).ravel()*0.1
#B = np.eye(num_groups)*0.85
#B = B + np.random.random(size=[num_groups,num_groups])*0.1

B = np.eye(num_groups)*0.8
B = B + np.ones([num_groups,num_groups])*0.2-np.eye(num_groups)*0.2

#---------------------------- Setup Model -----------------------------------#
raw_model = model.create_model(data_matrix, num_people, num_groups, alpha, B)
#model_instance = pymc.Model(raw_model)

#---------------------------- Call MAP to initialize MCMC -------------------#
#pymc.MAP(model_instance).fit(method='fmin_powell')
print '---------- Finished Running MAP to Set MCMC Initial Values ----------'
#---------------------------- Run MCMC --------------------------------------#
print '--------------------------- Starting MCMC ---------------------------'
M = pymc.MCMC(raw_model)
M.sample(100000,50000, thin=5, verbose=0)

