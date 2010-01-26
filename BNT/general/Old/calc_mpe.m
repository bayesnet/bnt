function [mpe, ll] = calc_mpe(engine, evidence, break_ties)
% CALC_MPE Computes the most probable explanation of the evidence
% [mpe, ll] = calc_mpe_given_inf_engine(engine, evidence, break_ties)
%
% INPUT
% engine must support max-propagation
% evidence{i} is the observed value of node i, or [] if hidden
% break_ties is optional. If 1, we will force ties to be broken consistently
%  by calling enter_evidence N times.
%
% OUTPUT
% mpe{i} is the most likely value of node i (cell array!)
% ll is the log-likelihood of the globally best assignment
%
% This currently only works when all hidden nodes are discrete

if nargin < 3, break_ties = 0; end


[engine, ll] = enter_evidence(engine, evidence, 'maximize', 1);

observed = ~isemptycell(evidence);

if 0 % fgraphs don't support bnet_from_engine
onodes = find(observed);
bnet = bnet_from_engine(engine);
pot_type = determine_pot_type(bnet, onodes);
assert(pot_type == 'd');
end

scalar = 1;
evidence = evidence(:); % hack to handle unrolled DBNs
N = length(evidence);
mpe = cell(1,N);
for i=1:N
  m = marginal_nodes(engine, i);
  % observed nodes are all set to 1 inside the inference engine, so we must undo this
  if observed(i)
    mpe{i} = evidence{i};
  else
    mpe{i} = argmax(m.T);
    % Bug fix by Ron Zohar, 8/15/01
    % If there are ties, we must break them as follows (see Jensen96, p106)
    if break_ties
      evidence{i} = mpe{i};                             
      [engine, ll] = enter_evidence(engine, evidence, 'maximize', 1);  
    end
  end
  if length(mpe{i}) > 1, scalar = 0; end
end

if nargout >= 2
  bnet = bnet_from_engine(engine);
  ll = log_lik_complete(bnet, mpe(:));
end
if 0 % scalar
  mpe = cell2num(mpe);
end
