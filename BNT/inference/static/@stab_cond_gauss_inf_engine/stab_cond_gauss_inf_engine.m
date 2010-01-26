function engine = stab_cond_gauss_inf_engine(bnet)
% STAB_COND_GAUSS_INF_ENGINE Junction tree using stable CG potentials
% engine = cond_gauss_inf_engine(bnet)
% 
% This class was written by Shan Huang (shan.huang@intel.com) 2001
% and fixed by Rainer Deventer deventer@informatik.uni-erlangen.de March 2003
N = length(bnet.dag);
clusters = {};
root = N;
stages = { 1:N };
onodes = [];
engine = init_fields;
engine.evidence = [];
engine = class(engine, 'stab_cond_gauss_inf_engine', inf_engine(bnet));

ns = bnet.node_sizes(:);
ns(onodes) = 1; % observed nodes have only 1 possible value

%[engine.jtree, dummy, engine.cliques, B, w, elim_order, moral_edges, fill_in_edges, strong] = ...
%    dag_to_jtree(bnet, onodes, stages, clusters);


partial_order = determine_elim_constraints(bnet, onodes);
strong = ~isempty(partial_order);
stages = {};
clusters = {};
[engine.jtree, dummy_root, engine.cliques, B, w, elim_order] = 
    graph_to_jtree(moralize(bnet.dag), ns, partial_order, stages, clusters);

    
engine.cliques_bitv = B;
engine.clique_weight = w;
C = length(engine.cliques);
engine.clpot = cell(1,C);

% A node can be a member of many cliques, but is assigned to exactly one, to avoid
% double-counting its CPD. We assign node i to clique c if c is the "lightest" clique that
% contains i's family, so it can accomodate its CPD.

engine.clq_ass_to_node = zeros(1, N);
num_cliques = length(engine.cliques);
for i=1:N
  clqs_containing_family = find(all(B(:,family(bnet.dag, i)), 2)); % all selected columns must be 1
  c = clqs_containing_family(argmin(w(clqs_containing_family)));  
  engine.clq_ass_to_node(i) = c; 
end

% Compute the separators between connected cliques.
[is,js] = find(engine.jtree > 0);
engine.separator = cell(num_cliques, num_cliques);
for k=1:length(is)
  i = is(k); j = js(k);
  engine.separator{i,j} = find(B(i,:) & B(j,:)); % intersect(cliques{i}, cliques{j});
end
%keyboard;
engine.seppot = cell(C,C);

pot_type = 'scg';
check_for_cd_arcs([], bnet.cnodes, bnet.dag);

% Make the jtree rooted, so there is a fixed message passing order.
if strong
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Start the search for the strong root at the clique with the  %
  % highest number.                                              %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  root = length(engine.cliques);
  root_found = 0;
  
  while ((~root_found) & (root >= 1))
      root_found = test_strong_root(engine.jtree,engine.cliques,bnet.dnodes,root);
      if ~root_found
          root = root - 1;
      end
  end
  assert(root > 0)
  engine.root = root;
  % the last clique is guaranteed to be a strong root
  %engine.root = length(engine.cliques);
else
  % jtree_dbn_inf_engine requires the root to contain the interface.
  % This may conflict with the strong root requirement! *********** BUG *************
  engine.root = clq_containing_nodes(engine, root);
  if engine.root <= 0
    error(['no clique contains ' num2str(root)]);
  end
end  

[engine.jtree, engine.preorder, engine.postorder] = mk_rooted_tree(engine.jtree, engine.root);

% Evaluate CPDs with evidence, and convert to potentials  
pot = cell(1, N);
inited = zeros(1, C);
clpot = cell(1, C);
evidence = cell(1, N);
for n=1:N
  fam = family(bnet.dag, n);
  e = bnet.equiv_class(n);
  %pot{n} = CPD_to_scgpot(bnet.CPD{e}, fam, ns, bnet.cnodes, evidence);
  pot{n} = convert_to_pot(bnet.CPD{e}, pot_type, fam(:), evidence);
  cindex = engine.clq_ass_to_node(n);
  if inited(cindex)
      clpot{cindex} = direct_combine_pots(pot{n}, clpot{cindex});
  else
      clpot{cindex} = pot{n};
      inited(cindex) = 1;
  end
end

for i=1:C
    if inited(i) == 0
        clpot{i} = scgpot([], [], [], []);
    end
end

seppot = cell(C, C);
% separators are is not need to initialize

% collect to root (node to parents)
% Unlike the HUGIN architecture the complements are stored in the cliques during COLLECT 
% and the separators are not playing a specific role during this process
for n=engine.postorder(1:end-1)
  for p=parents(engine.jtree, n)
    if ~isempty(engine.separator{p,n})
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % The empty case might happen for unlinked nodes, i.e. the DAG is not %
      % a single tree, but a forest                                           %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      [margpot, comppot] = complement_pot(clpot{n}, engine.separator{p,n});
      clpot{n} = comppot;
      clpot{p} = combine_pots(clpot{p}, margpot);
    end
  end
end

% distribute message from root
% We have not to store the weak clique marginals and keep the original complement potentials. 
% This is a minor variation of HUGIN architecture.
temppot = clpot;
for n=engine.preorder
  for c=children(engine.jtree, n)
    seppot{n,c} = marginalize_pot(temppot{n}, engine.separator{n,c});
    temppot{c} = direct_combine_pots(temppot{c}, seppot{n,c});
  end
end

engine.clpot = clpot;
engine.seppot = seppot;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init_fields()                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function engine = init_fields()

engine.evidence = [];
engine.jtree = [];
engine.cliques = [];
engine.cliques_bitv = [];
engine.clique_weight = [];
engine.preorder = [];
engine.postorder = [];
engine.root = []; 
engine.clq_ass_to_node = [];
engine.separator = [];
engine.clpot =[];
engine.seppot = [];












