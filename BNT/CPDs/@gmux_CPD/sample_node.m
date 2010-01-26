function y = sample_node(CPD, pev)
% SAMPLE_NODE Draw a random sample from P(Xi | x(pi_i), theta_i)  (gmux)
% y = sample_node(CPD, parent_evidence)
%
% parent_ev{i} is the value of the i'th parent

dpval = pev{CPD.dps};
x = pev{CPD.cps(dpval)};
y = gsamp(CPD.mean(:,dpval) + CPD.weights(:,:,dpval)*x(:), CPD.cov(:,:,dpval), 1);
y = y(:);
