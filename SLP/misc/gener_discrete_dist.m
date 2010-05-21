function beta = gener_discrete_dist(N, alpha)
% beta = gener_beta(N, masse)
%
%   using the logistic function to create a random vector
%   of discrete probability distribution on N values
%   if masse~=1 then sum(beta)==masse and this is not a probability
% 
% francois.olivier.c.h@gmail.com

if nargin==1, alpha=1; end
gamma = rand(1,N);
beta  = alpha*N*(exp(gamma)/sum(exp(gamma)));
