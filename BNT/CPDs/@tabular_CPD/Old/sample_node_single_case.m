function y = sample_node(CPD, pev)
% SAMPLE_NODE Draw a random sample from P(Xi | x(pi_i), theta_i)  (tabular)
% y = sample_node(CPD, pev)
%
% pev{i} is the value of the i'th parent (if any)

%assert(~any(isemptycell(pev)));

%CPT = CPD_to_CPT(CPD);
%sz = mysize(CPT);
sz = CPD.sizes; 
nparents = length(sz)-1;
if nparents > 0
  pvals = cat(1, pev{:});
end
switch nparents
 case 0, T = CPD.CPT;
 case 1, T = CPD.CPT(pvals(1), :);
 case 2, T = CPD.CPT(pvals(1), pvals(2), :);
 case 3, T = CPD.CPT(pvals(1), pvals(2), pvals(3), :);
 case 4, T = CPD.CPT(pvals(1), pvals(2), pvals(3), pvals(4), :);
 otherwise,
  psz = sz(1:end-1);
  ssz = sz(end);
  i = subv2ind(psz, pvals(:)');
  T = reshape(CPD.CPT, [prod(psz) ssz]);
  T = T(i,:);
end

if sz(end)==2
  r = rand(1,1);
  if r > T(1)
    y = 2;
  else
    y = 1;
  end
else
  y = sample_discrete(T);
end
