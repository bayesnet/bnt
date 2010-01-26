function y = sample_node(CPD, pev)
% SAMPLE_NODE Draw a random sample from P(Y|pa(y), theta)  (root)
% Y = SAMPLE_NODE(CPD, PEV)
%
% pev{i} is the evidence on the i'th parent.
% Since a root has no parents, we ignore pev,
% and return the value the root was clamped to when it was created.

y = CPD.val;
