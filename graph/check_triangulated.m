function [triangulated, order] = check_triangulated(G)
% CHECK_TRIANGULATED Return 1 if G is a triangulated (chordal) graph, 0 otherwise.
% [triangulated, order] = check_triangulated(G)
% 
% A numbering alpha is perfect if Nbrs(alpha(i)) intersect {alpha(1)...alpha(i-1)} is complete.
% A graph is triangulated iff it has a perfect numbering.
% The Maximum Cardinality Search algorithm will create such a perfect numbering if possible.
% See Golumbic, "Algorithmic Graph Theory and Perfect Graphs", Cambridge Univ. Press, 1985, p85.
% or Castillo, Gutierrez and Hadi, "Expert systems and probabilistic network models", Springer 1997, p134.


G = setdiag(G, 1);
n = length(G);
order = zeros(1,n);
triangulated = 1;
numbered = [1];
order(1) = 1;
for i=2:n
  U = mysetdiff(1:n, numbered); % unnumbered nodes
  score = zeros(1, length(U));
  for ui=1:length(U)
    u = U(ui);
    score(ui) = length(myintersect(neighbors(G, u), numbered));
  end
  u = U(argmax(score));
  numbered = [numbered u];
  order(i) = u;
  nns = myintersect(neighbors(G,u), order(1:i-1)); % numbered neighbors
  if ~isequal(G(nns,nns), ones(length(nns))) % ~complete(G(nns,nns))
    triangulated = 0;
    break;
  end
end

