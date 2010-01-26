function engine = gibbs_sampling_inf_engine(bnet, varargin)
% GIBBS_SAMPLING_INF_ENGINE
%
% engine = gibbs_sampling_inf_engine(bnet, ...) 
%
% Optional parameters [default in brackets]
% 'burnin' - How long before you start using the samples [100].
% 'gap' - how often you use the samples in the estimate [1].
% 'T' - number of samples [1000]
%   i.e, number of node flips (so, for
%   example if there are 10 nodes in the bnet, and T is 1000, each
%   node will get flipped 100 times (assuming a deterministic schedule)) 
%   The total running time is proportional to burnin + T*gap.
%
% 'order' - if the sampling schedule is deterministic, use this
% parameter to specify the order in which nodes are sampled.
% Order is allowed to include multiple copies of nodes, which is
% useful if you want to, say, focus sampling on particular nodes.
% Default is to use a deterministic schedule that goes through the
% nodes in order.
%
% 'sampling_dist' - when using a stochastic sampling method, at
% each step the node to sample is chosen according to this
% distribution (may be unnormalized)
% 
% The sampling_dist and order parameters shouldn't both be used,
% and this will cause an assert.
%
%
% Written by "Bhaskara Marthi" <bhaskara@cs.berkeley.edu> Feb 02.


engine.burnin = 100;
engine.gap = 1;
engine.T = 1000; 
use_default_order = 1;
engine.deterministic = 1;
engine.order = {};
engine.sampling_dist = {};

if nargin >= 2
  args = varargin;
  nargs = length(args);
  for i = 1:2:nargs
    switch args{i}
     case 'burnin'
      engine.burnin = args{i+1};
     case 'gap'
      engine.gap = args{i+1};
     case 'T'
      engine.T = args{i+1};
     case 'order'
      assert (use_default_order);
      use_default_order = 0;
      engine.order = args{i+1};
     case 'sampling_dist'
      assert (use_default_order);
      use_default_order = 0;
      engine.deterministic = 0;
      engine.sampling_dist = args{i+1};
     otherwise
      error(['unrecognized parameter to gibbs_sampling_inf_engine']);
    end
  end
end

engine.slice_size = size(bnet.dag, 2);
if (use_default_order)
  engine.order = 1:engine.slice_size;
end
engine.hnodes = [];
engine.onodes = [];
engine.evidence = [];
engine.state = [];
engine.marginal_counts = {};

% Precompute the strides for each CPT
engine.strides = compute_strides(bnet);

% Precompute graphical information
engine.families = compute_families(bnet);
engine.children = compute_children(bnet);

% For convenience, store the CPTs as tables rather than objects
engine.CPT = get_cpts(bnet);

engine = class(engine, 'gibbs_sampling_inf_engine', inf_engine(bnet));

















