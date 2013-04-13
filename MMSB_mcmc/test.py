import pymc
import numpy as np


data_matrix=np.loadtxt("../data/Y_alpha0.1_K5_N20.txt",delimiter=',')
num_people = 20
num_groups = 5
alpha = np.ones(num_groups).ravel()*0.1
B = np.eye(num_groups)*0.85
B = B + np.random.random(size=[num_groups,num_groups])*0.1


#---------------------------- Data Transform -----------------------------#
indices = np.indices(dimensions = np.shape(data_matrix))
people_indices = indices[0].reshape(num_people*num_people)
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

completed_pi_list = [pymc.CompletedDirichlet('completed_pi_%d' % i, dist) for i, dist in enumerate(pi_list)] 

# Indicator variables of whether the pth person is in a group or not
# DIMENSIONS: 1 x (num_people^2) for each list, where each element is Kx1
# DOMAIN : {0,1}, only one element of vector is 1, all else 0
# DISTRIBUTION: Categorical (using Multinomial with 1 observation)
z_pTq_matrix = np.empty([num_people,num_people], dtype=object)
z_pFq_matrix = np.empty([num_people,num_people], dtype=object)
for p_person in range(num_people):
    for q_person in range(num_people):
        z_pTq_matrix[p_person,q_person] = pymc.Multinomial('z_%dT%d_vector' % (p_person,q_person), n=1, p=pi_list[p_person], trace=False)
        z_pFq_matrix[p_person,q_person] = pymc.Multinomial('z_%dF%d_vector' % (p_person,q_person), n=1, p=pi_list[q_person], trace=False)

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
