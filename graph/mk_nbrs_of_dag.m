function [Gs, op, nodes] = mk_nbrs_of_dag(G0)
% MK_NBRS_OF_DAG Make all DAGs that differ from G0 by a single edge deletion, addition or reversal
% [Gs, op, nodes] = mk_nbrs_of_dag(G0)
%
% Gs{i} is the i'th neighbor.
% op{i} = 'add', 'del', or 'rev' is the operation used to create the i'th neighbor.
% nodes(i,1:2) are the head and tail of the operated-on arc.

Gs = {};
op = {};
nodes = [];

[I,J] = find(G0);
nnbrs = 1;
% all single edge deletions
for e=1:length(I)
  i = I(e); j = J(e);
  G = G0;
  G(i,j) = 0;
  Gs{nnbrs} = G;
  op{nnbrs} = 'del';
  nodes(nnbrs, :) = [i j];
  nnbrs = nnbrs + 1;
end

% all single edge reversals
for e=1:length(I)
  i = I(e); j = J(e);
  G = G0;
  G(i,j) = 0;
  G(j,i) = 1;
  if acyclic(G)
    Gs{nnbrs} = G;
    op{nnbrs} = 'rev';
    nodes(nnbrs, :) = [i j];
    nnbrs = nnbrs + 1;
  end
end

[I,J] = find(~G0);
% all single edge additions
for e=1:length(I)
  i = I(e); j = J(e);
  if i ~= j % don't add self arcs
    G = G0;
    G(i,j) = 1;
    if G(j,i)==0 % don't add i->j if j->i exists already
      if acyclic(G)
	Gs{nnbrs} = G;
	op{nnbrs} = 'add';
	nodes(nnbrs, :) = [i j];
	nnbrs = nnbrs + 1;
      end
    end
  end
end








