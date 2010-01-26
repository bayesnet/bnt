function y = sample_node(CPD, pev)
% SAMPLE_NODE Draw a random sample from P(Xi | x(pi_i), theta_i)  (gaussian)
% y = sample_node(CPD, parent_evidence)
%
% pev{i} is the value of the i'th parent (if there are any parents)
% y is the sampled value (a scalar or vector)

if length(CPD.dps)==0
  i = 1;
else
  dpvals = cat(1, pev{CPD.dps});
  i = subv2ind(CPD.sizes(CPD.dps), dpvals(:)');
end

if length(CPD.cps) == 0 
  y = gsamp(CPD.mean(:,i), CPD.cov(:,:,i), 1);
else
  pev = pev(:);
  x = cat(1, pev{CPD.cps});
  y = gsamp(CPD.mean(:,i) + CPD.weights(:,:,i)*x(:), CPD.cov(:,:,i), 1);
end
y = y(:);
