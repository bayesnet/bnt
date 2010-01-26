function [scenarios, log_probs] = enumerate_scenarios(bnet, evidence)
% ENUMERATE_SCENARIOS Enumerate all assignments, and return the prob. of the non-zeros ones
% function [scenarios, log_probs] = enumerate_scenarios(bnet, evidence)

assert(isempty(bnet.cnodes));
n = length(bnet.dag);
observed = ~isemptycell(evidence);
vals = cat(1,evidence{observed});
vals = vals(:)';
ns = bnet.node_sizes;

log_probs = [];
scenarios = [];
for i=1:prod(ns)
  inst = ind2subv(ns, i); % i'th instantiation
  if isempty(vals) | inst(observed) == vals % agrees with evidence
    ll = log_lik_complete(bnet, num2cell(inst(:)));
    log_probs = [log_probs ll];
    scenarios = [scenarios(:)' inst];
  end
end
