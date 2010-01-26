function [engine, loglik] = enter_evidence(engine, evidence, varargin)
% ENTER_EVIDENCE Add the specified evidence to the network (hmm)
% [engine, loglik] = enter_evidence(engine, evidence, ...)
%
% evidence{i,t} = [] if if X(i,t) is hidden, and otherwise contains its observed value (scalar or column vector)
%
% The following optional arguments can be specified in the form of name/value pairs:
% [default value in brackets]
%
% maximize - if 1, does max-product (not yet supported), else sum-product [0]
% filter   - if 1, does filtering, else smoothing [0]
% oneslice - 1 means only compute marginals on nodes within a single slice [0]
%
% e.g., engine = enter_evidence(engine, ev, 'maximize', 1)

maximize = 0;
filter = 0;
oneslice = 0;

% parse optional params
args = varargin;
nargs = length(args);
if nargs > 0
  for i=1:2:nargs
    switch args{i},
     case 'maximize', maximize = args{i+1}; 
     case 'filter',  filter = args{i+1}; 
     case 'oneslice', oneslice = args{i+1};
     otherwise,  
      error(['invalid argument name ' args{i}]);       
    end
  end
end

[ss T] = size(evidence);
engine.maximize = maximize;
engine.evidence = evidence;
bnet = bnet_from_engine(engine);
engine.node_sizes = repmat(bnet.node_sizes_slice(:), [1 T]);

obs_bitv = ~isemptycell(evidence(:));
bitv = reshape(obs_bitv, ss, T);
for t=1:T
  onodes = find(bitv(:,t));
  if ~isequal(onodes, bnet.observed(:))
    error(['dbn was created assuming observed nodes per slice were '...
	   num2str(bnet.observed(:)')  ' but the evidence in slice ' num2str(t) ...
	   ' has observed nodes ' num2str(onodes(:)')]);
  end
end

obslik = mk_hmm_obs_lik_matrix(engine, evidence);

%[alpha, beta, gamma, loglik, xi] = fwdback(engine.startprob, engine.transprob, obslik, ...
[alpha, beta, gamma, loglik, xi] = fwdback_twoslice(engine, engine.startprob,...
                                                    engine.transprob, obslik, ...
                                                    'maximize', maximize, 'fwd_only', filter, ...
                                                    'compute_xi', ~oneslice);

engine.one_slice_marginal = gamma; % gamma(:,t) for t=1:T
if ~oneslice
  Q = size(gamma,1);
  engine.two_slice_marginal = reshape(xi, [Q*Q T-1]); % xi(:,t) for t=1:T-1
end
