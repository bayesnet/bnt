function [engine, loglik] = enter_evidence(engine, evidence)
% ENTER_EVIDENCE Add the specified evidence to the network (gibbs_sampling_inf_engine)
% [engine, loglik] = enter_evidence(engine, evidence)
%
% evidence{i} = [] if if X(i) is hidden, and otherwise contains its observed value 
%
% loglik is not computed... we just return a 0 value

bnet = bnet_from_engine(engine);

engine.hnodes = find(isemptycell(evidence));
engine.onodes = mysetdiff(1:length(evidence), engine.hnodes);

engine.evidence = zeros(engine.slice_size, 1);

% Reset all counts since they are no longer valid
engine.marginal_counts = {};
%engine.state = sample_bnet (bnet, 1, 0);
engine.state = cell2num(sample_bnet(bnet));

% For speed, we use a normal (not cell) array.  We're making use of
% the current restriction to discrete nodes.
for i = engine.onodes
    engine.evidence(i) = evidence{i};
end

loglik = 0;


