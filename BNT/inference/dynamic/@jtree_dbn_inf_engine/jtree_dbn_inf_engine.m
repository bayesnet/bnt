function engine = jtree_dbn_inf_engine(bnet, varargin)
% JTREE_DBN_INF_ENGINE Junction tree inference algorithm for DBNs.
% engine = jtree_inf_engine(bnet, ...)
%
% The following optional arguments can be specified in the form of name/value pairs:
% [default value in brackets]
%
% clusters - specifies variables that must be grouped in the 1.5 slice DBN
% maximize - 1 means max-product, 0 means sum-product [0]
%
% e.g., engine = jtree_dbn_inf_engine(dbn, 'clusters', {[1 2]});
%
% This uses all of slice t-1 plus the backwards interface of slice t.
% By contrast, jtree_2TBN_inf_engine in the online directory uses
% the forwards interface of slice t-1 plus all of slice t.
% See my thesis for details.

ss = length(bnet.intra);

engine.maximize = 0;
clusters = {};

args = varargin;
for i=1:2:length(args)
  switch args{i},
   case 'clusters', clusters = args{i+1};
   case 'maximize', engine.maximize = args{i+1};
   otherwise, error(['unrecognized argument ' args{i}])
  end
end


engine.evidence = [];
engine.node_sizes = [];

[int, engine.persist, engine.transient] = compute_interface_nodes(bnet.intra, bnet.inter);
engine.interface = int;
engine.nonint = mysetdiff(1:ss, int);

onodes = bnet.observed;

if 0
  % Create a 2 slice jtree
  % We force there to be cliques containing the in and out interfaces for slices t and t+1.
  obs_nodes = [onodes(:) onodes(:)+ss];
  engine.jtree_engine = jtree_inf_engine(bnet, 'observed', obs_nodes(:), ...
					 'clusters', {int, int+ss}, 'root', int+ss);
else
  % Create a "1.5 slice" jtree, containing slice 1 and the interface nodes of slice 2
  % To keep the node numbering the same, we simply disconnect the non-interface nodes
  % from slice 2, and set their size to 1.
  % We do this to speed things up, and so that the likelihood is computed correctly - we do not need to do
  % this if we just want to compute marginals. 
  intra15 = bnet.intra;
  for i=engine.nonint(:)'
    intra15(i,:) = 0;
    intra15(:,i) = 0;
  end
  dag15 = [bnet.intra bnet.inter;
	 zeros(ss)    intra15];
  ns = bnet.node_sizes(:);
  ns(engine.nonint+ss) = 1; % disconnected nodes get size 1
  obs_nodes = [onodes(:) onodes(:)+ss];
  bnet15 = mk_bnet(dag15, ns, 'discrete', bnet.dnodes, 'equiv_class', bnet.equiv_class(:), ...
		   'observed', obs_nodes(:));

  %bnet15 = mk_dbn(intra15, bnet.inter, bnet.node_sizes_slice, bnet.dnodes_slice, ...
  %		  bnet.equiv_class(:,1), bnet.equiv_class(:,2), bnet.intra);
  % with the dbn, we can't independently control the sizes of slice 2 nodes
  
  if 1
    % use unconstrained elimination,
    % but force there to be a clique containing both interfaces
    clusters(end+1:end+2) = {int, int+ss};
    engine.jtree_engine = jtree_inf_engine(bnet15, 'clusters', clusters, 'root', int+ss);
  else
    % Use constrained elimination - this induces a clique that contain the 2nd interface,
    % but not the first.
    % Hence we throw in the first interface as an extra.
    stages = {1:ss, [1:ss]+ss};
    clusters(end+1:end+2) = {int, int+ss};
    engine.jtree_engine = jtree_inf_engine(bnet15, 'clusters', clusters, ...
					   'stages', stages, 'root', int+ss);
  end
end

engine.in_clq = clq_containing_nodes(engine.jtree_engine, int);
engine.out_clq = clq_containing_nodes(engine.jtree_engine, int+ss);
engine.jtree_struct = struct(engine.jtree_engine); % violate object privacy


% Also create an engine just for slice 1
bnet1 = mk_bnet(bnet.intra1, bnet.node_sizes_slice, 'discrete', myintersect(bnet.dnodes,1:ss), ...
		'equiv_class', bnet.equiv_class(:,1), 'observed', onodes);
for i=1:max(bnet1.equiv_class)
  bnet1.CPD{i} = bnet.CPD{i};
end

engine.jtree_engine1 = jtree_inf_engine(bnet1, 'clusters', {int}, 'root', int);

engine.in_clq1 = clq_containing_nodes(engine.jtree_engine1, int);
engine.jtree_struct1 = struct(engine.jtree_engine1); % violate object privacy

% stuff needed by marginal_nodes
engine.clpot = [];
engine.T = [];

engine = class(engine, 'jtree_dbn_inf_engine', inf_engine(bnet));

