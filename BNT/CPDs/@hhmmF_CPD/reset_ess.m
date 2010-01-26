function CPD = reset_ess(CPD)
% RESET_ESS Reset the Expected Sufficient Statistics of a hhmm F node.
% CPD = reset_ess(CPD)

CPD.sub_CPD_term = reset_ess(CPD.sub_CPD_term);   
