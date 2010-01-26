function CPD = update_ess(CPD, fmarginal, evidence, ns, cnodes, hidden_bitv)
% UPDATE_ESS Update the Expected Sufficient Statistics of a hhmmF node.
% function CPD = update_ess(CPD, fmarginal, evidence, ns, cnodes, hidden_bitv)
%
% We assume the F nodes are always hidden

% Figure out the node numbers associated with each parent
dom = fmarginal.domain;
%Fself = dom(end); 
%Fbelow = dom(CPD.Fbelow_ndx);
Qself = dom(CPD.Qself_ndx);
Qps = dom(CPD.Qps_ndx);

Qsz = CPD.Qsz;
Qpsz = CPD.Qpsz;

if all(hidden_bitv(Qps)) % we assume all are hidden or all are observed
  k_ndx = 1:Qpsz;
  eff_Qpsz = Qpsz;
else
  k_ndx = subv2ind(Qpsz, cat(1, evidence{Qps}));
  eff_Qpsz = 1;
end

if hidden_bitv(Qself)
  j_ndx = 1:Qsz;
  eff_Qsz = Qsz;
else
  j_ndx = evidence{Qself};
  eff_Qsz = 1;
end

% Fmarginal(Qps, Q, Fbelow, F)
fmarg = myreshape(fmarginal.T, [eff_Qpsz eff_Qsz  2 2]);

counts = zeros(Qpsz, Qsz, 2);
%counts(k_ndx, j_ndx, :) = sum(fmarginal.T(:, :, :, :), 3); % sum over Fbelow
counts(k_ndx, j_ndx, :) = fmarg(:, :, 2, :); % Fbelow = 2

CPD.sub_CPD_term = update_ess_simple(CPD.sub_CPD_term, counts);
