function CPD = learn_params(CPD,fam,data,ns,cnodes)
%LEARN_PARAMS Summary of this function goes here
%   Detailed explanation goes here

ncases = size(data, 2);
CPD = reset_ess(CPD);
% make a fully observed joint distribution over the family
fmarginal.domain = fam;
fmarginal.T = 1;
fmarginal.mu = [];
fmarginal.Sigma = [];
if ~iscell(data)
  cases = num2cell(data);
else
  cases = data;
end
hidden_bitv = zeros(1, max(fam));
for m=1:ncases
  % specify (as a bit vector) which elements in the family domain are hidden
  hidden_bitv = zeros(1, max(fmarginal.domain));
  ev = cases(:,m);
  hidden_bitv(find(isempty(ev)))=1;
  CPD = update_ess(CPD, fmarginal, ev, ns, cnodes, hidden_bitv);
end
CPD = maximize_params(CPD);
end

