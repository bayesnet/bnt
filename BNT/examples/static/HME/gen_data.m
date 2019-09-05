function [data, ndata1, ndata2, targets]=gen_data(ndata, seed)
% Generate data from three classes in 2d
% Setting 'seed' for reproducible results
% OUTPUT
% data : data set
% ndata1, ndata2: separator

if nargin<1,
   error('Missing data size');
end

input_dim = 2;
num_classes = 3;

if nargin==2,
   % Fix seeds for reproducible results
   randn('state', seed);
   rand('state', seed);
end

% Generate mixture of three Gaussians in two dimensional space
data = randn(ndata, input_dim);
targets = zeros(ndata, 3);

% Priors for the clusters
prior(1) = 0.4;
prior(2) = 0.3;
prior(3) = 0.3;

% Cluster centres
c = [2.0, 2.0; 0.0, 0.0; 1, -1];

ndata1 = round(prior(1)*ndata);
ndata2 = round((prior(1) + prior(2))*ndata);
% Put first cluster at (2, 2)
data(1:ndata1, 1) = data(1:ndata1, 1) * 0.5 + c(1,1);
data(1:ndata1, 2) = data(1:ndata1, 2) * 0.5 + c(1,2);
targets(1:ndata1, 1) = 1;

% Leave second cluster at (0,0)
data((ndata1 + 1):ndata2, :) = data((ndata1 + 1):ndata2, :);
targets((ndata1+1):ndata2, 2) = 1;

data((ndata2+1):ndata, 1) = data((ndata2+1):ndata,1) *0.6 + c(3, 1);
data((ndata2+1):ndata, 2) = data((ndata2+1):ndata,2) *0.6 + c(3, 2);
targets((ndata2+1):ndata, 3) = 1;

if 0
  ndata = 1;
  data = x;
  targets = [1 0 0];
end
