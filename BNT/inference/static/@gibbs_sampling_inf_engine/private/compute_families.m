function families = compute_families(bnet)
% COMPUTE_FAMILIES 
% precomputes the families of nodes in a bnet
%
% The return value is a cell array for now

ss = size(bnet.dag, 1);
families = cell(ss, 1);
for i = 1:ss
  families{i} = family(bnet.dag, i);
end

