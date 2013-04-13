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
    #indices = np.indices(dimesions = np.shape(data_matrix))
    #people_indices = indices[0].reshape(num_people*num_people)
    data_vector = data_matrix.reshape(num_people*num_people).T

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
    # DIMENSIONS: 1 x (num_people * num_groups)
    # SUPPORT: (0,1], Elements of each vector should sum to 1 for each person
    # DISTRIBUTION: Dirichlet(alpha)
    pi_list = np.empty(num_people, dtype=object)
    for person in range(num_people):
        person_pi = pymc.Dirichlet('pi_%i' % person, theta=alpha_vector)
        pi_list[person] = person_pi


    # Indicator variables of whether the pth person is in a group or not
    # DIMENSIONS: 1 x (num_people^2) for each list, where each element is Kx1
    # DOMAIN : {0,1}, only one element of vector is 1, all else 0
    # DISTRIBUTION: Categorical (using Multinomial with 1 observation)
    z_pTq_matrix = np.empty([num_people,num_people], dtype=object)
    z_pFq_matrix = np.empty([num_people,num_people], dtype=object)
    for p_person in range(num_people):
        for q_person in range(num_people):
            z_pTq_matrix[p_person,q_person] = pymc.Multinomial('z_%dT%d_vector' % (p_person,q_person), n=1, p=pi_list[p_person])
            z_pFq_matrix[p_person,q_person] = pymc.Multinomial('z_%dF%d_vector' % (p_person,q_person), n=1, p=pi_list[q_person])

    #---------------------------- Data Level ---------------------------------#
    # Combination of Priors to build the scalar parameter for y~Bernoulli
    @pymc.deterministic
    def bernoulli_parameters(z_pTq=z_pTq_matrix, z_pFq=z_pFq_matrix, B=B_matrix):
        """
        Takes in the two z_lists of Categorical Stochastic Objects
        Take their values (using Deterministic class)
        Dot Product with z'Bz
        """
        bernoulli_parameters = np.empty([num_people, num_people], dtype=object)
        for p in range(num_people):
            for q in range(num_people):
                bernoulli_parameters[p,q] = np.dot(np.dot(z_pTq[p,q], B), z_pFq[p,q])
        return bernoulli_parameters.reshape(1,num_people*num_people)

    # Observed response when person p is asked whether q is "connected"
    # Reshaped such that each person is asked sequentiall about all others, then
    # next person's vector, etc.
    #
    # Includes information about both p and q's group membership
    # y = Bern(z_p2p * B * z_q2p
    # y in {0,1}
    # DIMENSIONS: 1 x (num_people * num_people)
    y_vector = pymc.Bernoulli('y_vector', p=bernoulli_parameters, value=data_vector, observed=True)

    #---------------------------- Return all MCMC Objects --------------------#
    return locals()
