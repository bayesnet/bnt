function y = sample_node(CPD, pvals)
% SAMPLE_NODE Draw a random sample from P(Xi | x(pi_i), theta_i)  (discrete)
% y = sample_node(CPD, parent_evidence)
%
% parent_evidence{i} is the value of the i'th parent

if 0
n = length(pvals)+1;
dom = 1:n;
evidence = cell(1,n);
evidence(1:n-1) = pvals;
T = convert_to_table(CPD, dom, evidence);
y = sample_discrete(T);
end


CPT = CPD_to_CPT(CPD);
sz = mysize(CPT);
nparents = length(sz)-1;
switch nparents
 case 0, T = CPT;
 case 1, T = CPT(pvals{1}, :);
 case 2, T = CPT(pvals{1}, pvals{2}, :);
 case 3, T = CPT(pvals{1}, pvals{2}, pvals{3}, :);
 case 4, T = CPT(pvals{1}, pvals{2}, pvals{3}, pvals{4}, :);
 otherwise,
  pvals = cat(1, pvals{:});
  psz = sz(1:end-1);
  ssz = sz(end);
  i = subv2ind(psz, pvals(:)');
  T = reshape(CPT, [prod(psz) ssz]);
  T = T(i,:);
end
y = sample_discrete(T);
