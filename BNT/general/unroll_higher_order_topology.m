function M = unroll_higher_order_topology(intra, inter, T, intra1)
% UNROLL_DBN_TOPOLOGY Make the block diagonal adjacency matrix for a DBN consisting of T slices
% M = unroll_dbn_topology(intra, inter, T, intra1)
%
% intra is the connectivity within a slice, inter between two slices.
% M will have intra along the diagonal, and inter one above the diagonal.
% intra1 is an optional argumnet, in case the intra is different for the first slice.

if nargin < 4 
    intra1 = intra; 
end;


ss = length(intra); % slice size
M = sparse(ss*T, ss*T);
[rows,columns,order] = size(inter);
for t1 = 1:T
  b = 1 + (t1 - 1)*ss : t1*ss;
  if t1 == 1
      M(b,b) = intra1;
  else
      M(b,b) = intra;
  end
  for t2 = 1:order
    if t1 + t2 <= T
      M(b,b+t2*ss) = inter(:,:,t2);
    end
  end
end

