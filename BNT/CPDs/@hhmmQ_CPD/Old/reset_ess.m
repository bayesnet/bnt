function CPD = reset_ess(CPD)
% RESET_ESS Reset the Expected Sufficient Statistics of a hhmm Q node.
% CPD = reset_ess(CPD)

if ~isempty(CPD.sub_CPD_start)
  CPD.sub_CPD_start = reset_ess(CPD.sub_CPD_start);
end
CPD.sub_CPD_trans = reset_ess(CPD.sub_CPD_trans);   
