function engine = belprop_inf_engine(fg, max_iter, momentum, tol, maximize)

if nargin < 2, max_iter = length(fg.G); end
if nargin < 3, momentum = 0; end
if nargin < 4, tol = 1e-3; end
if nargin < 5, maximize = 0; end

engine.fgraph = fg;
engine.max_iter = max_iter;
engine.momentum = momentum;
engine.tol = tol;
engine.maximize = maximize;

% store results computed by enter_evidence here
ndoms = length(fg.doms);
nvars = length(fg.vars);
engine.marginal_domains = cell(1, ndoms);

% to compute the marginal on each variable, we need to know which domain to marginalize
% so we represent each domain as a bit vector, and compute its (pre-evidence) weight
engine.dom_weight = [];

% engine.dom_bitv = sparse(ndoms, nvars);
% ns = fg.node_sizes;
% for i=1:ndoms
%   engine.dom_bitv(i, fg.doms{i}) = 1;
%   engine.dom_weight(i) = prod(ns(fg.doms{i}));
% end


engine = class(engine, 'belprop_inf_engine');
