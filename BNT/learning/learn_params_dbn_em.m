function [bnet, LL, engine] = learn_params_dbn_em(engine, evidence, varargin)
% LEARN_PARAMS_DBN Set the parameters in a DBN to their ML/MAP values using batch EM.
% [bnet, LLtrace, engine] = learn_params_dbn_em(engine, data, ...)
%
% data{l}{i,t} = value of node i in slice t of time-series l, or [] if hidden.
%   Suppose you have L time series, each of length T, in an O*T*L array D, 
%   where O is the num of observed scalar nodes, and N is the total num nodes per slice.
%   Then you can create data as follows, where onodes is the index of the observable nodes:
%      data = cell(1,L);
%      for l=1:L
%        data{l} = cell(N, T);
%        data{l}(onodes,:) = num2cell(D(:,:,l));
%      end
% Of course it is possible for different sets of nodes to be observed in
% each slice/ sequence, and for each sequence to be a different length.
%
% LLtrace is the learning curve: the vector of log-likelihood scores at each iteration.
%
% Optional arguments [default]
%
% max_iter - specifies the maximum number of iterations [100]
% thresh - specifies the thresold for stopping EM [1e-3]
%   We stop when |f(t) - f(t-1)| / avg < threshold,
%   where avg = (|f(t)| + |f(t-1)|)/2 and f is log lik.
% verbose - display loglik at each iteration [1]
% anneal - 1 means do deterministic annealing (only for entropic priors) [0]
% anneal_rate - geometric cooling rate [0.8]
% init_temp - initial annealing temperature [10]
% final_temp - final annealing temperature [1e-3]
%

max_iter = 100;
thresh = 1e-3;
anneal = 0;
anneal_rate = 0.8;
init_temp = 10;
final_temp = 1e-3;
verbose = 1;

for i=1:2:length(varargin)
  switch varargin{i}
   case 'max_iter', max_iter = varargin{i+1}; 
   case 'thresh', thresh = varargin{i+1}; 
   case 'anneal', anneal = varargin{i+1}; 
   case 'anneal_rate', anneal_rate = varargin{i+1}; 
   case 'init_temp', init_temp = varargin{i+1}; 
   case 'final_temp', final_temp = varargin{i+1}; 
   otherwise, error(['unrecognized argument' varargin{i}])
  end
end

% take 1 EM step at each temperature value, then when temp=0, run to convergence
% When using an entropic prior, Z = 1-T, so 
% T=2 => Z=-1 (max entropy)
% T=1 => Z=0 (max likelihood)
% T=0 => Z=1 (min entropy / max structure)
num_iter = 1;
LL = [];
if anneal
  temperature = init_temp;
  while temperature > final_temp
    [engine, loglik, logpost] = EM_step(engine, evidence, temperature);
    if verbose
      fprintf('EM iteration %d, loglik = %8.4f, logpost = %8.4f, temp=%8.4f\n', ...
	      num_iter, loglik, logpost, temperature);
    end
    num_iter = num_iter + 1;
    LL = [LL loglik];
    temperature = temperature * anneal_rate;
  end
  temperature = 0;
  previous_loglik = loglik;
  previous_logpost = logpost;
else
  temperature = 0;
  previous_loglik = -inf;
  previous_logpost = -inf;
end

converged = 0;
while ~converged & (num_iter <= max_iter)
  [engine, loglik, logpost] = EM_step(engine, evidence, temperature);
  if verbose
    %fprintf('EM iteration %d, loglik = %8.4f, logpost = %8.4f\n', ...
    %	    num_iter, loglik, logpost);
    fprintf('EM iteration %d, loglik = %8.4f\n', num_iter, loglik);
  end
  num_iter = num_iter + 1;
  [converged, decreased] = em_converged(loglik, previous_loglik, thresh);
  %[converged, decreased] = em_converged(logpost, previous_logpost, thresh);
  previous_loglik = loglik;
  previous_logpost = logpost;
  LL = [LL loglik];
end

bnet = bnet_from_engine(engine);

%%%%%%%%%

function [engine, loglik, logpost] = EM_step(engine, cases, temp)

bnet = bnet_from_engine(engine); % engine contains the old params that are used for the E step
ss = length(bnet.intra);
CPDs = bnet.CPD; % these are the new params that get maximized
num_CPDs = length(CPDs);

% log P(theta|D) = (log P(D|theta) + log P(theta)) - log(P(D))
% where log P(D|theta) = sum_cases log P(case|theta)
% and log P(theta) = sum_CPDs log P(CPD) - only count once even if tied!
% logpost = log P(theta,D) (un-normalized)
% This should be negative, and increase at every step.

adjustable = zeros(1,num_CPDs);
logprior = zeros(1, num_CPDs);
for e=1:num_CPDs
  adjustable(e) = adjustable_CPD(CPDs{e});
end
adj = find(adjustable);

for e=adj(:)'
  logprior(e) = log_prior(CPDs{e});
  CPDs{e} = reset_ess(CPDs{e});
end

loglik = 0;
for l=1:length(cases)
  evidence = cases{l};
  if ~iscell(evidence)
    error('training data must be a cell array of cell arrays')
  end
  [engine, ll] = enter_evidence(engine, evidence);
  assert(~isnan(ll))
  loglik = loglik + ll;
  T = size(evidence, 2);
  
  % We unroll ns etc because in update_ess, we refer to nodes by their unrolled number
  % so that they extract evidence from the right place.
  % (The CPD should really store its own version of ns and cnodes...)
  ns = repmat(bnet.node_sizes_slice(:), [1 T]);
  cnodes = unroll_set(bnet.cnodes_slice, ss, T);
 
  %hidden_bitv = repmat(bnet.hidden_bitv(1:ss), [1 T]);
  hidden_bitv = zeros(ss, T);
  hidden_bitv(isemptycell(evidence))=1;
  % hidden_bitv(i) = 1 means node i is hidden.
  % We pass this in, rather than using isemptycell(evidence(dom)), because
  % isemptycell is very slow.
  
  t = 1;
  for i=1:ss
    e = bnet.equiv_class(i,1);
    if adjustable(e)
      fmarg = marginal_family(engine, i, t);
      CPDs{e} = update_ess(CPDs{e}, fmarg, evidence, ns(:), cnodes(:), hidden_bitv(:)); 
    end
  end
  
  for i=1:ss   
    e = bnet.equiv_class(i,2);
    if adjustable(e)
      for t=2:T
	fmarg = marginal_family(engine, i, t);
	CPDs{e} = update_ess(CPDs{e}, fmarg, evidence, ns(:), cnodes(:), hidden_bitv(:));
      end
    end
  end
end

logpost = loglik + sum(logprior(:));

for e=adj(:)'
  CPDs{e} = maximize_params(CPDs{e}, temp);
end

engine = update_engine(engine, CPDs);




