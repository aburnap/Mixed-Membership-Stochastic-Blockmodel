import numpy as np
import pymc
###############################################################################
#			Bayesian Model Structure Code for Mixed Membership Block Model
#			Authors: Alex Burnap, Efren Cruz, Xin Rong, Brian Segal
#			Date: April 1, 2013
###############################################################################

def create_model(data_matrix, num_people, num_groups, alpha, B):
    """
    Function that takes in a set of data and returns a dictionary of MCMC object
    of each random variable.  In addition, the hyperparameters for each random
    variable will be imbued in the objects.

    Overcommented for readibility's sake
    """
    #---------------------------- Data Transform -----------------------------#
    indices = np.indices(dimesions = np.shape(data_matrix))
    people_indices = indices[0].reshape(num_people*num_people)
    data_vector = data_matrix.reshape(num_people*num_groups)

    #---------------------------- Hyperparameters ----------------------------#
    # Average probability distribution of being in groups
    # DIMENSIONS: 1 x num_groups
    # SUPPORT: (0,inf)
    # DISTRIBUTION: None
    alpha_vector = alpha

    # Matrix of inter-group correlations
    # DIMENSIONS: num_groups x num_groups
    # SUPPORT: [0,1]
    # DISTRIBUTION: None
    B_matrix = B

    #---------------------------- Prior Parameters ---------------------------#
    # Actual group membership probabilities for each person
    # DIMENSIONS: num_people x num_groups
    # SUPPORT: (0,1], Elements of each vector should sum to 1 for each person
    # DISTRIBUTION: Dirichlet(alpha)
    pi_matrix = 
    
    # Indicator variable of whether the pth person is in a group or not
    # DIMENSIONS: 1 x (num_groups * num_people^2)
    # DOMAIN : {0,1}, only one element of vector is 1, all else 0
    # DISTRIBUTION: Categorical
    z_p2q_vector = pymc.Categorical('z_p2q_vector'
    z_q2p_vector
    
    #---------------------------- Data Level ---------------------------------#
    # Combination of Priors to build the scalar parameter for y~Bernoulli
    @pymc.deterministic
    def how_much_do_they_interact(z_p2q=z_p2q_vector, z_q2p=z_q2p_vector, B=B_vector):
        temp = np.dot(np.transpose(z_p2q), B)
        result = np.dot(temp, z_q2p)
        return np.array(result)
    
    def dot_product(x
    # Observed response when person p is asked whether q is "connected"
    # Reshaped such that each person is asked sequentiall about all others, then
    # next person's vector, etc.
    #
    # Includes information about both p and q's group membership
    # y = Bern(z_p2p * B * z_q2p
    # y in {0,1}
    # DIMENSIONS: 1 x (num_people * num_people)
    y_vector = pymc.Bernoulli('y_vector', PARAMETERS_HERE, value=data_vector, observed=True)
