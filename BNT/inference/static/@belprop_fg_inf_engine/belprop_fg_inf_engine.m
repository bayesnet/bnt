function engine = belprop_fg_inf_engine(fg, varargin) 
% BELPROP_FG_INF_ENGINE Make a belief propagation inference engine for factor graphs
% engine = belprop_fg_inf_engine(factor_graph, ...)
%
% The following optional arguments can be specified in the form of name/value pairs:
% [default in brackets]
% e.g., engine = belprop_inf_engine(fg, 'tol', 1e-2, 'max_iter', 10)
%
% max_iter - max. num. iterations [ 2*num_nodes ]
% momentum - weight assigned to old message in convex combination (useful for damping oscillations) [0]
% tol - tolerance used to assess convergence [1e-3]
% maximize - 1 means use max-product, 0 means use sum-product [0]
%
% This uses potential objects, like belprop_inf_engine, and hence is quite slow.

engine = init_fields;
engine = class(engine, 'belprop_fg_inf_engine');

% set params to default values
N = length(fg.G);
engine.max_iter = 2*N;
engine.momentum = 0;
engine.tol = 1e-3;
engine.maximize = 0;

% parse optional arguments
engine = set_params(engine, varargin);

engine.fgraph = fg;

% store results computed by enter_evidence here
engine.marginal_nodes = cell(1, fg.nvars);
engine.evidence = [];


%%%%%%%%%%%%

function engine = init_fields()

engine.fgraph = [];
engine.max_iter = [];
engine.momentum = [];
engine.tol = [];
engine.maximize = [];
engine.marginal_nodes = [];
engine.evidence = [];
engine.niter = [];
