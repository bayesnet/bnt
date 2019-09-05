function [M, moral_edges] = moralize(G)
% MORALIZE Ensure that for every child, all its parents are married, and drop directionality of edges.
% [M, moral_edges] = moralize(G)

M = G;
n = length(M);
for i=1:n
  fam = family(G,i);
  M(fam,fam)=1;
end
M = setdiag(M,0);
moral_edges = sparse(triu(max(0,M-G),1));
