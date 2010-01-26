function [engine, loglik] = enter_evidence(engine, evidence, varargin)
% ENTER_EVIDENCE Add the specified evidence to the network (var_elim)
% [engine, loglik] = enter_evidence(engine, evidence, ...)
%
% evidence{i} = [] if if X(i) is hidden, and otherwise contains its observed value (scalar or column vector)

% we could pre-process the evidence here, to prevent repeated work, but we don't.
engine.evidence = evidence;

if nargout == 2
  [m, loglik] = marginal_nodes(engine, [1]);
end
