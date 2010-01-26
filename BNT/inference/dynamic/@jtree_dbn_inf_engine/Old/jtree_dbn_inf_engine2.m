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


% Create a 2 slice jtree
% We force there to be cliques containing the in and out interfaces for slices t and t+1.
obs_nodes = [onodes(:) onodes(:)+ss];
engine.jtree_engine = jtree_inf_engine(bnet, 'observed', obs_nodes(:), ...
					 'clusters', {int, int+ss}, 'root', int+ss);

engine.in_clq = clq_containing_nodes(engine.jtree_engine, int);
engine.out_clq = clq_containing_nodes(engine.jtree_engine, int+ss);
engine.jtree_struct = struct(engine.jtree_engine); % violate object privacy



% Also create an engine just for slice 1
bnet1 = mk_bnet(bnet.intra1, bnet.node_sizes_slice, bnet.dnodes, bnet.equiv_class(:,1));
for i=1:max(bnet1.equiv_class)
  bnet1.CPD{i} = bnet.CPD{i};
end

engine.jtree_engine1 = jtree_inf_engine(bnet1, 'observed', onodes, 'clusters', {int}, ...
					'root', int);

engine.in_clq1 = clq_containing_nodes(engine.jtree_engine1, int);
engine.jtree_struct1 = struct(engine.jtree_engine1); % violate object privacy




% stuff needed by marginal_nodes
engine.clpot = [];
engine.T = [];
engine.maximize = [];

engine = class(engine, 'jtree_dbn_inf_engine', inf_engine(bnet));

