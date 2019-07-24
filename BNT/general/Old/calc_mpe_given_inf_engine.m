function [mpe, prob] = calc_mpe_given_inf_engine(engine, evidence)
% CALC_MPE_GIVEN_ENGINE Computes the most probable explanation of the evidence
% [mpe, prob] = calc_mpe_given_inf_engine(engine, evidence)
%
% INPUT
% engine must support max-propagation
% evidence{i} is the obsevred value of node i, or [] if hidden
%
% OUTPUT
% mpe(i) is the most likely value of node i
% prob is the likelihood of the globally best assignment
%
% This currently only works when all nodes are discrete

[engine, ll] = enter_evidence(engine, evidence);

observed = ~isemptycell(evidence);
N = length(evidence);
mpe = zeros(1,N);
for i=1:N
  m = marginal_nodes(engine, i);
  % discrete observed nodes are all set to 1 inside the inference engine, so we must undo this
  if observed(i)
    mpe(i) = evidence{i};
  else
    mpe(i) = argmax(m.T);
  end
end

bnet = bnet_from_engine(engine);
ll = log_lik_complete(bnet, num2cell(mpe(:)));
prob = exp(ll);
