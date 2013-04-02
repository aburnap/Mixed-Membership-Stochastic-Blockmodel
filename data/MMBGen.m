function [ Y, Pi ] = MMBGen( alpha, Beta, N )
%MMBGen Generates data for Mixed Membership Stochastic Blockmodel
%simulation
%   alpha: 1-by-K vector, each component should be greater than 0
%   Beta: K-by-K matrix, each cell should be a scalar between [0,1]
%   N: number of nodes
%Returns:
%   Y: N-by-N observation matrix, each cell is a binary value
%   Pi: N-by-K ground-truth memberships for each node.
Y = zeros(N, N);
K = length(alpha);
Pi = zeros(N, K);
for i = 1:N
    X = gamrnd(alpha, ones(1,K));
    Pi(i,:) = X/sum(X);
end
for i = 1:N
    for j = 1:N
        zi = mnrnd(1, Pi(i,:))';
        zj = mnrnd(1, Pi(j,:))';
        Y(i,j) = binornd(1, zi'*Beta*zj);
    end
end
end