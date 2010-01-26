function M = unroll_dbn_topology(intra, inter, T, intra1)
% UNROLL_DBN_TOPOLOGY Make the block diagonal adjacency matrix for a DBN consisting of T slices
% M = unroll_dbn_topology(intra, inter, T, intra1)
%
% intra is the connectivity within a slice, inter between two slices.
% M will have intra along the diagonal, and inter one above the diagonal.
% intra1 is an optional argumnet, in case the intra is different for the first slice.

if nargin < 4, intra1 = intra; end

ss = length(intra); % slice size
M = sparse(ss*T, ss*T);

b = 1:ss;
M(b,b) = intra1;
M(b,b+ss) = inter;

for t=2:T-1
  b = (1:ss) + (t-1)*ss;
  M(b,b) = intra;
  M(b,b+ss) = inter;
end

t = T;
b = (1:ss) + (t-1)*ss;
M(b,b) = intra;

   
