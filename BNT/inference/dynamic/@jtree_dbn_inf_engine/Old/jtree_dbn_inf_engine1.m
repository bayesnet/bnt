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

if 1
  % Create a 2 slice jtree
  % We force there to be cliques containing the in and out interfaces for slices t and t+1.
  obs_nodes = [onodes(:) onodes(:)+ss];
  engine.jtree_engine = jtree_inf_engine(bnet, 'observed', obs_nodes(:), ...
					 'clusters', {int, int+ss}, 'root', int+ss);
else
  % Create a "1.5 slice" jtree, containing slice 1 and the interface nodes of slice 2
  % To keep the node numbering the same, we simply disconnect the non-interface nodes
  % from slice 2.
  intra15 = bnet.intra;
  for i=engine.nonint(:)'
    intra15(i,:) = 0;
    intra15(:,i) = 0;
  end
  bnet15 = mk_dbn(intra15, bnet.inter, bnet.node_sizes_slice, bnet.dnodes_slice, ...
		  bnet.equiv_class(:,1), bnet.equiv_class(:,2), bnet.intra);
  obs_nodes = [onodes(:) onodes(:)+ss];
  engine.jtree_engine = jtree_inf_engine(bnet15, 'observed', obs_nodes(:), ...
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

