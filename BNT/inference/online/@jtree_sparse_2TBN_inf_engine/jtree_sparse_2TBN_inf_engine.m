function engine = jtree_sparse_2TBN_inf_engine(bnet, varargin)
% JTREE_ONLINE_INF_ENGINE Online Junction tree inference algorithm for DBNs.
% engine = jtree_online_inf_engine(bnet, ...)
%
% The following optional arguments can be specified in the form of name/value pairs:
% [default value in brackets]
%
% clusters - specifies variables that must be grouped in the 1.5 slice DBN
% maximize - 1 means do max-product, 0 means sum-product [0]
%
% The same nodes must be observed in every slice.

ss = length(bnet.intra);
clusters = {};
engine.maximize = 0;

args = varargin;
nargs = length(args);
for i=1:2:length(args)
  switch args{i},
   case 'clusters', clusters = args{i+1};
   case 'maximize', engine.maximize = args{i+1};
   otherwise, error(['unrecognized argument ' args{i}])
  end
end

engine.evidence = [];
engine.node_sizes = [];

int = [];
% include nodes with any outgoing arcs
for u=1:ss
  if any(bnet.inter(u,:))
    int = [int u];
  end
end

engine.interface = int;
engine.nonint = mysetdiff(1:ss, int);

onodes = bnet.observed;

% Create a "1.5 slice" jtree, containing the interface nodes of slice 1
% and all the nodes of slice 2
% To keep the node numbering the same, we simply disconnect the non-interface nodes
% from slice 1, and set their size to 1.
% We do this to speed things up, and so that the likelihood is computed correctly - we do not need to do
% this if we just want to compute marginals (i.e., we can include nodes whose potentials will
% be left as all 1s).
intra15 = bnet.intra;
for i=engine.nonint(:)'
  intra15(:,i) = 0;
  intra15(i,:) = 0;
  assert(~any(bnet.inter(i,:)))
end
dag15 = [intra15      bnet.inter;
	 zeros(ss)    bnet.intra];
ns = bnet.node_sizes(:);
ns(engine.nonint) = 1; % disconnected nodes get size 1
obs_nodes = [onodes(:) onodes(:)+ss];
bnet15 = mk_bnet(dag15, ns, 'discrete', bnet.dnodes, 'equiv_class', bnet.equiv_class(:), ...
		 'observed', obs_nodes(:));

% use unconstrained elimination,
% but force there to be a clique containing both interfaces
clusters(end+1:end+2) = {int, int+ss};
%engine.jtree_engine = jtree_inf_engine(bnet15, 'clusters', clusters, 'root', int+ss);
engine.jtree_engine = jtree_sparse_inf_engine(bnet15, 'clusters', clusters, 'root', int+ss);
jtree_engine = struct(engine.jtree_engine); % violate object privacy

engine.in_clq = clq_containing_nodes(engine.jtree_engine, int);
engine.out_clq = clq_containing_nodes(engine.jtree_engine, int+ss);
engine.clq_ass_to_node = jtree_engine.clq_ass_to_node;
engine.root = jtree_engine.root_clq;

% Also create an engine just for slice 1
bnet1 = mk_bnet(bnet.intra1, bnet.node_sizes_slice, 'discrete', myintersect(bnet.dnodes,1:ss), ...
		'equiv_class', bnet.equiv_class(:,1), 'observed', onodes);
for i=1:max(bnet1.equiv_class)
  bnet1.CPD{i} = bnet.CPD{i};
end
%engine.jtree_engine1 = jtree_inf_engine(bnet1, 'clusters', {int}, 'root', int);
engine.jtree_engine1 = jtree_sparse_inf_engine(bnet1, 'clusters', {int}, 'root', int);
jtree_engine1 = struct(engine.jtree_engine1); % violate object privacy
engine.int_clq1 = clq_containing_nodes(engine.jtree_engine1, int);
engine.clq_ass_to_node1 = jtree_engine1.clq_ass_to_node;
engine.root1 = jtree_engine1.root_clq;

engine.observed = [onodes onodes+ss];
engine.observed1 = onodes;
engine.pot_type = determine_pot_type(bnet, onodes);
engine.slice_size = bnet.nnodes_per_slice;

engine = class(engine, 'jtree_sparse_2TBN_inf_engine', inf_engine(bnet));

