function [Gs, op, nodes] = mk_nbrs_of_digraph2(G0)
% MK_NBRS_OF_DIGRAPH Make all digraphs that differ from G0 by a single edge deletion, addition or reversal
% [Gs, op, nodes] = mk_nbrs_of_digraph(G0)
% op{i} = 'add', 'del', or 'rev' is the operation used to create the i'th neighbor.
% nodes(i,1:2) are the head and tail of the operated-on arc.

[I,J] = find(G0);
G0bar = setdiag(~G0, 0); % exclude self loops in graph complement
[Ibar,Jbar] = find(G0bar);
nnbrs = 2*length(I) + length(Ibar);
Gs = cell(1, nnbrs);
op = cell(1, nnbrs);
nodes = zeros(nnbrs, 2);

nbr = 1;
% all single edge deletions
for e=1:length(I)
  i = I(e); j = J(e);
  G = G0;
  G(i,j) = 0;
  Gs{nbr} = G;
  op{nbr} = 'del';
  nodes(nbr, :) = [i j];
  nbr = nbr + 1;
end

% all single edge reversals
for e=1:length(I)
  i = I(e); j = J(e);
  G = G0;
  G(i,j) = 0;
  G(j,i) = 1;
  Gs{nbr} = G;
  op{nbr} = 'rev';
  nodes(nbr, :) = [i j];
  nbr = nbr + 1;
end

[I,J] = find(~G0);
% all single edge additions
for e=1:length(I)
  i = I(e); j = J(e);
  G = G0;
  if i ~= j % don't add self loops
    G(i,j) = 1;
    Gs{nbr} = G;
    op{nbr} = 'add';
    nodes(nbr, :) = [i j];
    nbr = nbr + 1;
  end
end

assert(nnbrs == nbr-1);
