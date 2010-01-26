function engine = hmm_2TBN_inf_engine(bnet, varargin)
% HMM_2TBN_INF_ENGINE Inference engine for DBNs which uses the forwards-backwards algorithm.
% engine = hmm_2TBN_inf_engine(bnet, ...)
%
% The DBN is converted to an HMM with a single meganode, but the observed nodes remain factored.
% This can be faster than jtree if the num. hidden nodes is low, because of lower constant factors.
%
% All hidden nodes must be discrete.
% All observed nodes are assumed to be leaves.
% The parents of each observed leaf are assumed to be a subset of the hidden nodes within the same slice.
% The only exception is if bnet is an AR-HMM, where the parents are assumed to be self in the
% previous slice (continuous), plus all the discrete nodes in the current slice.


%% Optional arguments
%% ndx_type - 'B', 'D', or 'SD', used in marginal_family [ 'SD' ]

ndx_type = 'SD';
ss = bnet.nnodes_per_slice;

% parse optional params
args = varargin;
nargs = length(args);
if nargs > 0
  for i=1:2:nargs
    switch args{i},
     %case 'ndx_type', ndx_type = args{i+1};
     otherwise,  
      error(['invalid argument name ' args{i}]);       
    end
  end
end

% Stuff to do with speeding up marginal_family
%engine.ndx_type = ndx_type;

[int, engine.persist, engine.transient] = compute_interface_nodes(bnet.intra, bnet.inter);
engine.persist_bitv = zeros(1, ss);
engine.persist_bitv(engine.persist) = 1;


ns = bnet.node_sizes(:);
ns(bnet.observed) = 1;
ns(bnet.observed+ss) = 1;
engine.eff_node_sizes = ns;

% for n=1:ss
%   dom = 1:(2*ss); % domain of xi(:,:,1)
%   fam = family(bnet.dag, n+ss);
%   engine.marg_fam2_ndx_id(n) = add_ndx(dom, fam, ns, ndx_type);
 
%   dom = 1:ss; % domain of gamma(:,:,1)
%   fam = family(bnet.dag, n);
%   engine.marg_fam1_ndx_id(n) = add_ndx(dom, fam, ns, ndx_type);

%   engine.marg_singleton_ndx_id(n) = add_ndx(dom, n, ns, ndx_type);
% end

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
engine.maximize = [];
engine.evidence = [];
engine.node_sizes = [];

% avoid the need to do bnet_from_engine, which is slow
engine.slice_size = ss;
engine.parents = bnet.parents;

engine.bel = [];
engine = class(engine, 'hmm_2TBN_inf_engine', inf_engine(bnet));

