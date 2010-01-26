function families = compute_families_dbn(bnet)
% COMPUTE_FAMILIES 
% precomputes the families of nodes in a dbn
%
% The return value is a cell array for now

ss = size(bnet.intra, 1);
families = cell(ss, 2);
for i = 1:ss
  families{i, 1} = family(bnet.dag, i, 1);
  families{i, 2} = family(bnet.dag, i, 2);
end

