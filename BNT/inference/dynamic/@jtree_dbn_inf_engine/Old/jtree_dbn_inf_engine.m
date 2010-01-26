function engine = jtree_dbn_inf_engine(bnet, varargin)
% JTREE_DBN_INF_ENGINE Junction tree inference algorithm for DBNs.

ss = length(bnet.intra);

onodes = [];

if nargin >= 2
  args = varargin;
  nargs = length(args);
  for i=1:2:nargs
    switch args{i},
     case 'observed', onodes = args{i+1};
    end
  end
end

[int, engine.persist, engine.transient] = compute_interface_nodes(bnet.intra, bnet.inter);
%engine.interface = engine.persist; % WRONG!
engine.interface = int;
engine.nonint = mysetdiff(1:ss, int);

if 0
  % Create a 2 slice jtree
  % We force there to be cliques containing the in and out interfaces for slices t and t+1.
  obs_nodes = [onodes(:) onodes(:)+ss];
  engine.jtree_engine = jtree_inf_engine(bnet, 'observed', obs_nodes(:), ...
					 'clusters', {int, int+ss}, 'root', int+ss);
else
  % Create a "1.5 slice" jtree, containing slice 1 and the interface nodes of slice 2
  nodes15 = [1:ss int+ss];
  N = length(nodes15);
  dag15 = bnet.dag(nodes15, nodes15);
  ns15 = bnet.node_sizes(nodes15);
  eclass15 = bnet.equiv_class(nodes15);
  discrete_bitv = zeros(1,2*ss);
  discrete_bitv(bnet.dnodes) = 1;
  discrete15 = find(discrete_bitv(nodes15));
  bnet15 = mk_bnet(dag15, ns15, 'equiv_class', eclass15, 'discrete', discrete15);
  bnet15.CPD = bnet.CPD; % CPDs for non-interface nodes in slice 2 will not be used
  obs_bitv = zeros(1, 2*ss);
  obs_bitv([onodes onodes+ss]) = 1;
  obs_nodes15 = find(obs_bitv(nodes15));
  int_bitv = zeros(1,ss);
  int_bitv(int) = 1;
  engine.jtree_engine = jtree_inf_engine(bnet15, 'observed', obs_nodes15(:), ...
				     'clusters', {int, int+ss}, 'root', int+ss);
end

engine.in_clq = clq_containing_nodes(engine.jtree_engine, int);
engine.out_clq = clq_containing_nodes(engine.jtree_engine, int+ss);

engine.clq_ass_to_node = zeros(ss, 2);
for i=1:ss
  engine.clq_ass_to_node(i, 1) = clq_containing_nodes(engine.jtree_engine, i);
  engine.clq_ass_to_node(i, 2) = clq_containing_nodes(engine.jtree_engine, i+ss);
end

engine.jtree_struct = struct(engine.jtree_engine); % violate object privacy

% stuff needed by marginal_nodes
engine.clpot = [];
engine.maximize = [];
engine.T = [];

engine = class(engine, 'jtree_dbn_inf_engine', inf_engine(bnet));

