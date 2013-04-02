###############################################################################
#				Run Code for Mixed Membership Block Model
#				Authors: Alex Burnap, Efren Cruz, Xin Rong, Brian Segal
#				Date: April 1, 2013
###############################################################################


############################# Dependencies ####################################
import MMBM.model as model
import numpy as np
import matplotlib.pyplot as plt
from Brians_modus_operandi import jigga_fresh_nasty as cur_method

############################# Run Time Parameters #############################

# Probably going to try and use - flag conventions with __init__(self, *args, **kwargs)


############################# Load Data #######################################
X=np.loadtxt("../data/Pi_alpha0.1_K5_N20.txt")

