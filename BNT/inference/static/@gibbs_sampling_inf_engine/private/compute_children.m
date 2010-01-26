function c = compute_children(bnet)
% COMPUTE_CHILDREN
% precomputes the children of nodes in a bnet
%
% The return value is a cell array for now

ss = size(bnet.dag, 1);
c = cell(ss, 1);
for i = 1:ss
  c{i} = children(bnet.dag, i);
end

