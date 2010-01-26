function bnet = bayes_update_params(bnet, cases, clamped)
% BAYES_UPDATE_PARAMS Bayesian parameter updating given completely observed data
% bnet = bayes_update_params(bnet, cases, clamped)
%
% If there is a missing data, you must use EM.
% cases(i,m) is the value assigned to node i in case m (this can also be a cell array).
% clamped(i,m) = 1 if node i was set by intervention in case m (default: clamped = zeros).
% Clamped nodes are not updated.
% If there is a single case, clamped is a list of the clamped nodes, not a bit vector.


%if iscell(cases), usecell = 1; else usecell = 0; end

n = length(bnet.dag);
ncases = size(cases, 2);
if n ~= size(cases, 1)
  error('data must be of size nnodes * ncases');
end

if ncases == 1 % clamped is a list of nodes
  if nargin < 3, clamped = []; end
  clamp_set = clamped;
  clamped = zeros(n,1);
  clamped(clamp_set) = 1;
else % each row of clamped is a bit vector
  if nargin < 3, clamped = zeros(n,ncases); end
end

for i=1:n
  e = bnet.equiv_class(i);
  if adjustable_CPD(bnet.CPD{e})
    u = find(clamped(i,:)==0);
    ps = parents(bnet.dag, i);
    bnet.CPD{e} = bayes_update_params(bnet.CPD{e}, cases(i,u), cases(ps,u));
  end
end


