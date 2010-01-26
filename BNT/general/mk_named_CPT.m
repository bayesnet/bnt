function CPT2 = mk_named_CPT(family_names, names, dag, CPT1)
% MK_NAMED_CPT Permute the dimensions of a CPT so they agree with the internal numbering convention
% CPT2 = mk_named_CPT(family_names, names, dag, CPT1)
%
% This is best explained by example.
% Consider the following directed acyclic graph
%
%     C
%   /   \
%  R     S
%   \   /
%     W
%
% where all arcs point down.
% When we create the CPT for node W, we consider S as its first parent, and R as its
% second, and hence write
%
%      S R W
% CPT1(1,1,:) = [1.0 0.0];
% CPT1(2,1,:) = [0.2 0.8];  % P(W=1 | R=1, S=2) = 0.2
% CPT1(1,2,:) = [0.1 0.9]; 
% CPT1(2,2,:) = [0.01 0.99];
%
% However, when we create the dag using mk_adj_mat, the nodes get topologically sorted,
% and by chance, node R preceeds node S in this ordering.
% Hence we should have written
%
%      R S W
% CPT2(1,1,:) = [1.0 0.0];
% CPT2(2,1,:) = [0.1 0.9];
% CPT2(1,2,:) = [0.2 0.8]; % P(W=1 | R=1, S=2) = 0.2
% CPT2(2,2,:) = [0.01 0.99];
%
% Since we do not know the order of the nodes in advance, we can write
%   CPT2 = mk_named_CPT({'S', 'R', 'W'}, names, dag, CPT1)
% where 'S', 'R', 'W' are the order of the dimensions we assumed (the child node must be last in this list),
% and names{i} is the name of the i'th node.

n = length(family_names);
family_nums = zeros(1,n);
for i=1:n
  family_nums(i) = stringmatch(family_names{i}, names); % was strmatch
end

fam = family(dag, family_nums(end));
perm = zeros(1,n);
for i=1:n
  %  perm(i) = find(family_nums(i) == fam);
  perm(i) = find(fam(i) == family_nums);
end

CPT2 = permute(CPT1, perm);
