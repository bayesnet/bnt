function engine = jtree_limid_inf_engine(bnet)
% JTREE_LIMID_INF_ENGINE Make a junction tree engine for use by solve_limid
% engine = jtree_limid_inf_engine(bnet)
%
% This engine is designed to compute marginals on decision nodes


MG = moralize(bnet.dag);
% We do not remove the utility nodes, because that complicates the book-keeping.
% Leaving them in will not introduce any un-necessary triangulation arcs, because they are always leaves.
% Also, since utility nodes have size 1, they do not increase the size of the potentials.

ns = bnet.node_sizes;
elim_order = best_first_elim_order(MG, ns);
[MTG, engine.cliques]  = triangulate(MG, elim_order);
[engine.jtree, root, B, w] = cliques_to_jtree(engine.cliques, ns);

% A node can be a member of many cliques, but is assigned to exactly one, to avoid
% double-counting its CPD. We assign node i to clique c if c is the "lightest" clique that
% contains i's family, so it can accomodate its CPD.
N = length(bnet.dag);
engine.clq_ass_to_node = zeros(1, N);
for i=1:N
  clqs_containing_family = find(all(B(:,family(bnet.dag, i)), 2)); % all selected columns must be 1
  c = clqs_containing_family(argmin(w(clqs_containing_family)));  
  engine.clq_ass_to_node(i) = c; 
end


% Compute the separators between connected cliques.
[is,js] = find(engine.jtree > 0);
num_cliques = length(engine.cliques);
engine.separator = cell(num_cliques, num_cliques);
for k=1:length(is)
  i = is(k); j = js(k);
  engine.separator{i,j} = find(B(i,:) & B(j,:)); % intersect(cliques{i}, cliques{j});
end


% create |D| different rooted jtree's
engine.rooted_jtree = cell(1, N);
engine.preorder = cell(1, N);
engine.postorder = cell(1, N);
for d=bnet.decision_nodes(:)'
  root = engine.clq_ass_to_node(d);
  [engine.rooted_jtree{d}, engine.preorder{d}, engine.postorder{d}] = mk_rooted_tree(engine.jtree, root);
end

engine.exclude = [];
engine.evidence = [];

engine = class(engine, 'jtree_limid_inf_engine', inf_engine(bnet));
