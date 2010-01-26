function c = clq_containing_nodes(engine, nodes, fam)
% CLQ_CONTAINING_NODES Find the lightest clique (if any) that contains the set of nodes
% c = clq_containing_nodes(engine, nodes, family)
%
% If the optional 'family' argument is specified, it means nodes = family(nodes(end)).
% (This is useful since clq_ass_to_node is not accessible to outsiders.)
% Returns c=-1 if there is no such clique.

if nargin < 3, fam = 0; else fam = 1; end

if length(nodes)==1
  c = engine.clq_ass_to_node(nodes(1));
elseif fam
  c = engine.clq_ass_to_node(nodes(end));
else
  B = engine.cliques_bitv;
  w = engine.clique_weight;
  clqs = find(all(B(:,nodes), 2)); % all selected columns must be 1
  if isempty(clqs)
    c = -1;
  else
    c = clqs(argmin(w(clqs)));     
  end
end
