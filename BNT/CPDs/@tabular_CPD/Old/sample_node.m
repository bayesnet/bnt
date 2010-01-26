function y = sample_node(CPD, pev, nsamples)
% SAMPLE_NODE Draw a random sample from P(Xi | x(pi_i), theta_i)  (tabular)
% Y = SAMPLE_NODE(CPD, PEV, NSAMPLES)
%
% pev(i,m) is the value of the i'th parent in sample m (if there are any parents).
% y(m) is the m'th sampled value (a row vector).
% (If pev is a cell array, so is y.)
% nsamples defaults to 1.

if nargin < 3, nsamples = 1; end

%if nargin < 4, usecell = 0; end
if iscell(pev), usecell = 1; else usecell = 0; end

if nsamples == 1, pev = pev(:); end

sz = CPD.sizes; 
nparents = length(sz)-1;
if nparents==0
  y = sample_discrete(CPD.CPT, 1, nsamples);
  if usecell
    y = num2cell(y);
  end
  return;
end

sz = CPD.sizes; 
[nparents nsamples] = size(pev);

if usecell
  pvals = cell2num(pev)'; % each row is a case
else
  pvals = pev';
end

psz = sz(1:end-1);
ssz = sz(end);
ndx = subv2ind(psz, pvals);
T = reshape(CPD.CPT, [prod(psz) ssz]);
T2 = T(ndx,:); % each row is a distribution selected by the parents
C = cumsum(T2, 2); % sum across columns
R = rand(nsamples, 1);
y = ones(nsamples, 1);
for i=1:ssz-1
  y = y + (R > C(:,i));
end
y = y(:)';
if usecell
  y = num2cell(y);
end



