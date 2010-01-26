function y = sample_node(CPD, pvals)
% SAMPLE_NODE Draw a random sample from P(Xi | x(pi_i), theta_i)  (discrete)
% y = sample_node(CPD, parent_evidence)
%
% parent_evidence{i} is the value of the i'th parent

n = length(pvals)+1;
dom = 1:n;
%evidence = cell(1,n);
%evidence(1:n-1) = pvals(:)';
evidence = pvals;
evidence{end+1} = [];
T = convert_to_table(CPD, dom, evidence);
y = sample_discrete(T);
