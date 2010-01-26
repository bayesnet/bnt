function engine = jtree_2TBN_inf_engine(bnet, varargin)
% JTREE_ONLINE_INF_ENGINE Online Junction tree inference algorithm for DBNs.
% engine = jtree_online_inf_engine(bnet, ...)
%
% The following optional arguments can be specified in the form of name/value pairs:
% [default value in brackets]
%
% clusters - specifies variables that must be grouped in the 1.5 slice DBN
%
% The same nodes must be observed in every slice.
%
% This uses the forwards interface of slice t-1 plus all of slice t.
% By contrast, jtree_dbn uses all of slice t-1 plus the backwards interface of slice t.
% See my thesis for details.


clusters = {};

args = varargin;
nargs = length(args);
for i=1:2:length(args)
  switch args{i},
   case 'clusters', clusters = args{i+1};
   otherwise, error(['unrecognized argument ' args{i}])
  end
end

engine.maximize = 0;
engine.evidence = [];
engine.node_sizes = [];

int = compute_fwd_interface(bnet.intra, bnet.inter);
engine.interface = int;
ss = length(bnet.intra);
engine.nonint = mysetdiff(1:ss, int);
onodes = bnet.observed;

bnet15 = mk_slice_and_half_dbn(bnet, int);

% use unconstrained elimination,
% but force there to be a clique containing both interfaces
clusters(end+1:end+2) = {int, int+ss};
engine.jtree_engine = jtree_inf_engine(bnet15, 'clusters', clusters, 'root', int+ss);
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
engine.jtree_engine1 = jtree_inf_engine(bnet1, 'clusters', {int}, 'root', int);
jtree_engine1 = struct(engine.jtree_engine1); % violate object privacy
engine.int_clq1 = clq_containing_nodes(engine.jtree_engine1, int);
engine.clq_ass_to_node1 = jtree_engine1.clq_ass_to_node;
engine.root1 = jtree_engine1.root_clq;

engine.observed = [onodes onodes+ss];
engine.observed1 = onodes;
engine.pot_type = determine_pot_type(bnet, onodes);
engine.slice_size = bnet.nnodes_per_slice;

engine = class(engine, 'jtree_2TBN_inf_engine', inf_engine(bnet));

