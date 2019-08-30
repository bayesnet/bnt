function CPD = update_ess(CPD, fmarginal, evidence, ns, cnodes, hidden_bitv)
% UPDATE_ESS Update the Expected Sufficient Statistics of a Von Mises node
% function CPD = update_ess(CPD, fmarginal, evidence, ns, cnodes, hidden_bitv)

%if nargin < 6
%  hidden_bitv = zeros(1, max(fmarginal.domain));
%  hidden_bitv(find(isempty(evidence)))=1;
%end

dom = fmarginal.domain;
ps = dom(1:end-1);
cps = myintersect(ps, cnodes);

CPD.nsamples = CPD.nsamples + 1;            
[ss cpsz dpsz] = size(CPD.weights); % ss = self size
[ss dpsz] = size(CPD.mean);

% general (non-vectorized) case
fullm = add_evidence_to_vmarginal(fmarginal, evidence, ns, cnodes); % slow!

if dpsz == 1 % no discrete parents
  w = 1;
else
  w = fullm.T(:); 
end

CPD.Wsum = CPD.Wsum + w;
yi = (cpsz+1):(cpsz+ss);
for i=1:dpsz
  muY = fullm.mu(yi, i);
  CPD.WYsum(:,i) = CPD.WYsum(:,i) + w(i)*cos(muY); %SS for Von Mises
  CPD.WYYsum(:,:,i) = CPD.WYYsum(:,:,i) + w(i)*(sin(muY)); %SS for Von Mises
end                

