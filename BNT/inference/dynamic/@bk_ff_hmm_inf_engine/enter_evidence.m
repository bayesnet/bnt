function [engine, loglik] = enter_evidence(engine, evidence, varargin)
% ENTER_EVIDENCE Add the specified evidence to the network (bk_ff_hmm)
% [engine, loglik] = enter_evidence(engine, evidence, ...)
%
% evidence{i,t} = [] if if X(i,t) is hidden, and otherwise contains its observed value (scalar or column vector)
%
% The following optional arguments can be specified in the form of name/value pairs:
% [default value in brackets]
%
% maximize - if 1, does max-product (not yet supported), else sum-product [0]
% filter -   if 1, do filtering, else smoothing [0]
%
% e.g., engine = enter_evidence(engine, ev, 'maximize', 1)

maximize = 0;
filter = 0;

% parse optional params
args = varargin;
nargs = length(args);
if nargs > 0
  for i=1:2:nargs
    switch args{i},
     case 'maximize', maximize = args{i+1}; 
     case 'filter', filter = args{i+1}; 
     otherwise,  
      error(['invalid argument name ' args{i}]);       
    end
  end
end

assert(~maximize);

bnet = bnet_from_engine(engine);
ss = length(bnet.intra);
onodes = bnet.observed;
hnodes = mysetdiff(1:ss, onodes);
T = size(evidence, 2);
assertBNT(~any(isemptycell(evidence(onodes,:))));

obslik = mk_hmm_obs_lik_mat(bnet, onodes, evidence);

ns = bnet.node_sizes_slice;
ns(onodes) = 1;

[gamma, loglik, marginals, marginalsT] = bk_ff_fb(engine.prior, engine.transmat, obslik, filter, hnodes, ns);
  
for t=1:T
  for i=hnodes(:)'
    engine.marginals{i,t} = pot_to_marginal(marginalsT{i,t});
  end
  for i=onodes(:)'
    m.domain = i + (t-1)*ss;
    m.T = 1;
    engine.marginals{i,t} = m;
  end
end



