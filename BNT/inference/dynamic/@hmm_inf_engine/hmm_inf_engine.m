function engine = hmm_inf_engine(bnet, varargin)
% HMM_INF_ENGINE Inference engine for DBNs which uses the forwards-backwards algorithm.
% engine = hmm_inf_engine(bnet, ...)
%
% The following optional arguments can be specified in the form of name/value pairs:
% [default value in brackets]
%
% maximize - 1 means max-product, 0 means sum-product [0]
%
% The DBN is converted to an HMM with a single meganode, but the observed nodes remain factored.
% This can be faster than jtree if the num. hidden nodes is low, because of lower constant factors.
%
% All hidden nodes must be discrete.
% All observed nodes are assumed to be leaves, i.e., they cannot be parents of anything.
% The parents of each observed leaf are assumed to be a subset of the hidden nodes within the same slice.
% The only exception is if bnet is an AR-HMM, where the parents are assumed to be self in the
% previous slice (continuous), plus all the discrete nodes in the current slice.

ss = bnet.nnodes_per_slice;

engine.maximize = 0;
% parse optional params
args = varargin;
nargs = length(args);
if nargs > 0
  for i=1:2:nargs
    switch args{i},
     case 'maximize', engine.maximize = args{i+1};
     otherwise,  
      error(['invalid argument name ' args{i}]);       
    end
  end
end

% Stuff to do with speeding up marginal_family
[int, engine.persist, engine.transient] = compute_interface_nodes(bnet.intra, bnet.inter);
engine.persist_bitv = zeros(1, ss);
engine.persist_bitv(engine.persist) = 1;


ns = bnet.node_sizes(:);
ns(bnet.observed) = 1;
ns(bnet.observed+ss) = 1;
engine.eff_node_sizes = ns;

for o=bnet.observed(:)'
  %if bnet.equiv_class(o,1) ~= bnet.equiv_class(o,2)
  %  error(['observed node ' num2str(o) ' is not tied'])
  %end
  cs = children(bnet.dag, o);
  if ~isempty(cs)
    error(['observed node ' num2str(o) ' is not allowed children'])
  end
end

[engine.startprob, engine.transprob, engine.obsprob] = dbn_to_hmm(bnet);

% This is where we will store the results between enter_evidence and marginal_nodes
engine.one_slice_marginal = [];
engine.two_slice_marginal = [];

ss = length(bnet.intra);
engine.evidence = [];
engine.node_sizes = [];

% avoid the need to do bnet_from_engine, which is slow
engine.slice_size = ss;
engine.parents = bnet.parents;

engine = class(engine, 'hmm_inf_engine', inf_engine(bnet));

