function engine = jtree_mnet_inf_engine(model, varargin)
% JTREE_MNET_INF_ENGINE Junction tree inference engine for Markov nets
% engine = jtree_inf_engine(mnet, ...)
%

% set default params
N = length(mnet.graph);
root = N;

engine = init_fields;
engine = class(engine, 'jtree_mnet_inf_engine', inf_engine(bnet));

onodes = bnet.observed;
if is_mnet(bnet)
  MG = bnet.graph;
else
  error('should be a mnet')
end

%[engine.jtree, dummy, engine.cliques, B, w, elim_order, moral_edges, fill_in_edges, strong] = ...
%    dag_to_jtree(bnet, onodes, stages, clusters);

porder = determine_elim_constraints(bnet, onodes);
strong = ~isempty(porder);
ns = bnet.node_sizes(:);
ns(onodes) = 1; % observed nodes have only 1 possible value
[engine.jtree, root2, engine.cliques, B, w] = ...
    graph_to_jtree(MG, ns, porder, stages, clusters);

engine.cliques_bitv = B;
engine.clique_weight = w;
C = length(engine.cliques);
engine.clpot = cell(1,C);

% Compute the separators between connected cliques.
[is,js] = find(engine.jtree > 0);
engine.separator = cell(C,C);
for k=1:length(is)
  i = is(k); j = js(k);
  engine.separator{i,j} = find(B(i,:) & B(j,:)); % intersect(cliques{i}, cliques{j});
end

% A node can be a member of many cliques, but is assigned to exactly one, to avoid
% double-counting its CPD. We assign node i to clique c if c is the "lightest" clique that
% contains i's family, so it can accomodate its CPD.

engine.clq_ass_to_node = zeros(1, N);
for i=1:N
  %c = clq_containing_nodes(engine, family(bnet.dag, i));
  clqs_containing_family = find(all(B(:,family(bnet.dag, i)), 2)); % all selected columns must be 1
  c = clqs_containing_family(argmin(w(clqs_containing_family)));  
  engine.clq_ass_to_node(i) = c; 
end

% Make the jtree rooted, so there is a fixed message passing order.
if strong
  % the last clique is guaranteed to be a strong root
  engine.root_clq = length(engine.cliques);
else
  % jtree_dbn_inf_engine requires the root to contain the interface.
  % This may conflict with the strong root requirement! *********** BUG *************
  engine.root_clq = clq_containing_nodes(engine, root);
  if engine.root_clq <= 0
    error(['no clique contains ' num2str(root)]);
  end
end  

[engine.jtree, engine.preorder, engine.postorder] = mk_rooted_tree(engine.jtree, engine.root_clq);

% collect 
engine.postorder_parents = cell(1,length(engine.postorder));
for n=engine.postorder(:)'
  engine.postorder_parents{n} = parents(engine.jtree, n);
end
% distribute
engine.preorder_children = cell(1,length(engine.preorder));
for n=engine.preorder(:)'
  engine.preorder_children{n} = children(engine.jtree, n);
end

  

%%%%%%%%

function engine = init_fields()

engine.jtree = [];
engine.cliques = [];
engine.separator = [];
engine.cliques_bitv = [];
engine.clique_weight = [];
engine.clpot = [];
engine.clq_ass_to_node = [];
engine.root_clq = [];
engine.preorder = [];
engine.postorder = [];
engine.preorder_children = [];
engine.postorder_parents = [];
engine.maximize = [];
engine.evidence = [];

