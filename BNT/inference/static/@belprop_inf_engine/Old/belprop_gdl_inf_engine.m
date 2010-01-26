function engine = belprop_gdl_inf_engine(gdl, varargin) 
% BELPROP_GDL_INF_ENGINE Make a belief propagation inference engine for a GDL graph
% engine = belprop_gdl_inf_engine(gdl_graph, ...)
%
% If the GDL graph is a tree, this will give exact results.
%
% The following optional arguments can be specified in the form of name/value pairs:
% [default in brackets]
% e.g., engine = belprop_inf_engine(gdl, 'tol', 1e-2, 'max_iter', 10)
%
% protocol - 'tree' means send messages up then down the tree,
%            'parallel' means use synchronous updates ['parallel']
% max_iter - max. num. iterations [ 2*num_nodes ]
% momentum - weight assigned to old message in convex combination (useful for damping oscillations) [0]
% tol - tolerance used to assess convergence [1e-3]
% maximize - 1 means use max-product, 0 means use sum-product [0]


engine = init_fields;
engine = class(engine, 'belprop_gdl_inf_engine');

% set default params
N = length(gdl.G);
engine.protocol = 'parallel';
engine.max_iter = 2*N;
engine.momentum = 0;
engine.tol = 1e-3;
engine.maximize = 0;

engine = set_params(engine, varargin);

engine.gdl = gdl;

if strcmp(engine.protocol, 'tree')
  % Make a rooted tree, so there is a fixed message passing order.
  root = N;
  [engine.tree, engine.preorder, engine.postorder, height, cyclic] = mk_rooted_tree(gdl.G, root);
  assert(~cyclic);
end

% store results computed by enter_evidence here
ndoms = length(gdl.doms);
nvars = length(gdl.vars);
engine.marginal_domains = cell(1, ndoms);

% to compute the marginal on each variable, we need to know which domain to marginalize
% and we want to choose the lightest. We compute the weight once we have seen the evidence.
engine.dom_weight = [];
engine.evidence = [];


%%%%%%%%%

function engine = init_fields()

engine.protocol = [];
engine.gdl = [];
engine.max_iter = [];
engine.momentum = [];
engine.tol = [];
engine.maximize = [];
engine.marginal_domains = [];
engine.evidence = [];
engine.tree = [];
engine.preorder = [];
engine.postorder = [];
engine.dom_weight = [];
