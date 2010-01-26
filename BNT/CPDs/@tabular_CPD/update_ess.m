function CPD = update_ess(CPD, fmarginal, evidence, ns, cnodes, hidden_bitv)
% UPDATE_ESS Update the Expected Sufficient Statistics of a tabular node.
% function CPD = update_ess(CPD, fmarginal, evidence, ns, cnodes, hidden_bitv)

dom = fmarginal.domain;

if all(hidden_bitv(dom))
  CPD = update_ess_simple(CPD, fmarginal.T);
  %fullm = add_ev_to_dmarginal(fmarginal, evidence, ns);
  %assert(approxeq(fullm.T(:), fmarginal.T(:)))
else
  fullm = add_ev_to_dmarginal(fmarginal, evidence, ns);
  CPD = update_ess_simple(CPD, fullm.T);
end

