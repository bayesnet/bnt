function score = BIC_score_CPD(CPD, fam, data, ns, cnodes)
% BIC_score_CPD Compute the BIC score of a generic CPD
% score = BIC_score_CPD(CPD, fam, data, ns, cnodes)
%
% We assume this node has a maximize_params method

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
for m=1:ncases
  CPD = update_ess(CPD, fmarginal, cases(:,m), ns, cnodes);
end
CPD = maximize_params(CPD);
self = fam(end);
ps = fam(1:end-1);
L = log_prob_node(CPD, cases(self,:), cases(ps,:));
score = L - 0.5*CPD.nparams*log(ncases);
