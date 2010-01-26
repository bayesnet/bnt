function [engine, loglik] = enter_evidence(engine, evidence)
% ENTER_EVIDENCE Add the specified evidence to the network (enumerative_inf)
% [engine, loglik] = enter_evidence(engine, evidence)
%
% evidence{i} = [] if if X(i) is hidden, and otherwise contains its observed value (scalar or column vector)

engine.evidence = evidence;
if nargout == 2
  [m, loglik] = marginal_nodes(engine, []);
end
