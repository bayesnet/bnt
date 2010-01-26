function engine = bk_inf_engine(bnet, varargin)
% BK_INF_ENGINE Boyen-Koller approximate inference algorithm for DBNs.
%
% In the BK algorithm, the belief state is represented as a product of marginals,
% even though the factors may not be independent.
%
% engine = bk_inf_engine(bnet, ...)
%
% The following optional arguments can be specified in the form of name/value pairs:
% [default value in brackets]
% 
% clusters - if a cell array, clusters{i} specifies the terms in the i'th factor.
%          - 'exact' means create one cluster that contains all the nodes in a slice [exact]
%          - 'ff' means create one cluster per node (ff = fully factorised).
%
%
% For details, see
% - "Tractable Inference for Complex Stochastic Processes", X. Boyen and D. Koller, UAI 98.
% - "Approximate learning of dynamic models",  X. Boyen and D. Koller, NIPS 98.
% (The UAI98 paper discusses filtering and theory, and the NIPS98 paper discusses smoothing.)

ss = length(bnet.intra);
% set default params
clusters = 'exact';


if nargin >= 2
  args = varargin;
  nargs = length(args);
  for i=1:2:nargs
    switch args{i},
     case 'clusters',  clusters = args{i+1};
     otherwise, error(['unrecognized argument ' args{i}])
    end
  end
end

if strcmp(clusters, 'exact')
  %clusters = { compute_interface_nodes(bnet.intra, bnet.inter) };
  clusters = { 1:ss }; 
elseif strcmp(clusters, 'ff')
  clusters = num2cell(1:ss);
end


% We need to insert the prior on the clusters in slice 1,
% and extract the posterior on the clusters in slice 2.
C = length(clusters);
clusters2 = cell(1,2*C);
clusters2(1:C) = clusters;
for c=1:C
  clusters2{c+C} = clusters{c} + ss;
end

onodes = bnet.observed;
obs_nodes = [onodes(:) onodes(:)+ss];
engine.sub_engine = jtree_inf_engine(bnet, 'clusters', clusters2);

engine.clq_ass_to_cluster = zeros(C, 2);
for c=1:C
  engine.clq_ass_to_cluster(c,1) = clq_containing_nodes(engine.sub_engine, clusters{c});
  engine.clq_ass_to_cluster(c,2) = clq_containing_nodes(engine.sub_engine, clusters{c}+ss);
end
engine.clusters = clusters;

engine.clq_ass_to_node = zeros(ss, 2);
for i=1:ss
  engine.clq_ass_to_node(i, 1) = clq_containing_nodes(engine.sub_engine, i);
  engine.clq_ass_to_node(i, 2) = clq_containing_nodes(engine.sub_engine, i+ss);
end



% Also create an engine just for slice 1
bnet1 = mk_bnet(bnet.intra1, bnet.node_sizes_slice, 'discrete', myintersect(bnet.dnodes, 1:ss), ...
		'equiv_class', bnet.equiv_class(:,1), 'observed', onodes);
for i=1:max(bnet1.equiv_class)
  bnet1.CPD{i} = bnet.CPD{i};
end

engine.sub_engine1 = jtree_inf_engine(bnet1, 'clusters', clusters);

engine.clq_ass_to_cluster1 = zeros(1,C);
for c=1:C
  engine.clq_ass_to_cluster1(c) = clq_containing_nodes(engine.sub_engine1, clusters{c});
end

engine.clq_ass_to_node1 = zeros(1, ss);
for i=1:ss
  engine.clq_ass_to_node1(i) = clq_containing_nodes(engine.sub_engine1, i);
end

engine.clpot = []; % this is where we store the results between enter_evidence and marginal_nodes
engine.filter = [];
engine.maximize = [];
engine.T = [];

engine.bel = [];
engine.bel_clpot = [];
engine.slice1 = [];
%engine.pot_type = 'cg';
% hack for online inference so we can cope with hidden Gaussians and discrete
% it will not affect the pot type used in enter_evidence
engine.pot_type = determine_pot_type(bnet, onodes);

engine = class(engine, 'bk_inf_engine', inf_engine(bnet));

