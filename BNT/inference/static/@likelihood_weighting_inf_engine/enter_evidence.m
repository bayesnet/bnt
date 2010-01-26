function [engine, ll] = enter_evidence(engine, evidence, nsamples)
% ENTER_EVIDENCE Add the specified evidence to the network (likelihood_weighting)
% [engine, ll] = enter_evidence(engine, evidence, nsamples)
% evidence{i} = [] if if X(i) is hidden, and otherwise contains its observed value (scalar or column vector)
%
% If nsamples is not specified, the value specified when the engine was created will be used.
% ll (log-likelihood) is set to [].

ll = [];
if nargin < 3, nsamples = engine.nsamples; end

bnet = bnet_from_engine(engine);
N = length(bnet.dag);
samples = cell(nsamples, N);
weights = zeros(1, nsamples);

ns = bnet.node_sizes;
original_evidence = evidence;
observed = ~isemptycell(original_evidence);
for s=1:nsamples
  evidence = original_evidence(:); % must be a column vector
  w = 1;
  for i=1:N
    ps = parents(bnet.dag, i);
    e = bnet.equiv_class(i);
    if observed(i)
      p = exp(log_prob_node(bnet.CPD{e}, evidence(i), evidence(ps)));
      w = w * p;
    else
      x = sample_node(bnet.CPD{e}, evidence(ps));
      evidence{i} = x;
    end
  end
  samples(s,:) = evidence;
  weights(s) = w;
end                 

engine.samples = samples;
engine.weights = weights;
