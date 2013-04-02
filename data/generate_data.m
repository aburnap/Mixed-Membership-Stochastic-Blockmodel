% Generate data for MMB analysis
% Specify parameters
K = 5;
N = 20;
alpha_value = 0.1;

% Specifiy alpha parameter

alpha = ones(1,K)*alpha_value;

% Beta: inter-cluster connectivity matrix
% The diagonal values of Beta should be quite big (close to 1)
% The other elements should be close to 0
Beta = diag(unifrnd(ones(1,K)*0.7, ones(1,K)*0.9));
Beta = Beta + unifrnd(zeros(K,K), ones(K,K)*0.1);

[Y, Pi] = MMBGen( alpha, Beta, N );