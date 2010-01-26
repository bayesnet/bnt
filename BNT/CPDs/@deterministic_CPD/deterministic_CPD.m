function CPD = deterministic_CPD(bnet, self, fname, pfail)
% DETERMINISTIC_CPD Make a tabular CPD representing a (noisy) deterministic function
%
% CPD = deterministic_CPD(bnet, self, fname)
% This calls feval(fname, pvals) for each possible vector of parent values.
% e.g., suppose there are 2 ternary parents, then pvals = 
%  [1 1], [2 1], [3 1],   [1 2], [2 2], [3 2],   [1 3], [2 3], [3 3]
% If v = feval(fname, pvals(i)), then
%  CPD(x | parents=pvals(i)) = 1 if x==v, and = 0 if x<>v
% e.g., suppose X4 = X2 AND (NOT X3). Then
%    bnet.CPD{4} = deterministic_CPD(bnet, 4, inline('((x(1)-1) & ~(x(2)-1)) + 1'));  
% Note that x(1) refers pvals(1) = X2, and x(2) refers to pvals(2)=X3
% See also boolean_CPD.
%
% CPD = deterministic_CPD(bnet, self, fname, pfail)
% will put probability mass 1-pfail on f(parents), and distribute pfail over the other values.
% This is useful for simulating noisy deterministic functions.
% If pfail is omitted, it is set to 0.
%


if nargin==0
  % This occurs if we are trying to load an object from a file.
  CPD = tabular_CPD(bnet, self);
  return;
elseif isa(bnet, 'deterministic_CPD')
  % This might occur if we are copying an object.
  CPD = bnet;
  return;
end

if nargin < 4, pfail = 0; end

ps = parents(bnet.dag, self);
ns = bnet.node_sizes;
psizes = ns(ps);
self_size = ns(self);

psucc = 1-pfail;

CPT = zeros(prod(psizes), self_size);
pvals = zeros(1, length(ps));
for i=1:prod(psizes)
  pvals = ind2subv(psizes, i);
  x = feval(fname, pvals);
  %fprintf('%d ', [pvals x]); fprintf('\n');
  if psucc == 1
    CPT(i, x) = 1;
  else
    CPT(i, x) = psucc;
    rest = mysetdiff(1:self_size, x);
    CPT(i, rest) = pfail/length(rest);
  end
end
CPT = reshape(CPT, [psizes self_size]);  

CPD = tabular_CPD(bnet, self, 'CPT',CPT, 'clamped',1);


