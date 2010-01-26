function [B,D,mu] = extract_params_from_gbn(bnet)
% Extract all the local parameters of each Gaussian node, and collect them into global matrices.
% [B,D,mu] = extract_params_from_gbn(bnet)
%
% B(i,j) is a block matrix that contains the transposed weight matrix from node i to node j.
% D(i,i) is a block matrix that contains the noise covariance matrix for node i.
% mu(i) is a block vector that contains the shifted noise mean for node i.

% In Shachter's model, the mean of each node in the global gaussian is
% the same as the node's local unconditional mean.
% In Alag's model (which we use), the global mean gets shifted.


num_nodes = length(bnet.dag);
bs = bnet.node_sizes(:); % bs = block sizes
N = sum(bs); % num scalar nodes

B = zeros(N,N);
D = zeros(N,N);
mu = zeros(N,1);

for i=1:num_nodes % in topological order
  ps = parents(bnet.dag, i);
  e = bnet.equiv_class(i);
  %[m, Sigma, weights] = extract_params_from_CPD(bnet.CPD{e});
  s = struct(bnet.CPD{e}); % violate privacy of object
  m = s.mean; Sigma = s.cov; weights = s.weights;
  if length(ps) == 0
    mu(block(i,bs)) = m;
  else
    mu(block(i,bs)) = m + weights *  mu(block(ps,bs));
  end
  B(block(ps,bs), block(i,bs)) = weights';
  D(block(i,bs), block(i,bs)) = Sigma;
end



