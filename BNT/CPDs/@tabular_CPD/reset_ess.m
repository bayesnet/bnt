function CPD = reset_ess(CPD)
% RESET_ESS Reset the Expected Sufficient Statistics of a tabular node.
% CPD = reset_ess(CPD)

%CPD.counts = zeros(size(CPD.CPT));
CPD.counts = zeros(prod(size(CPD.CPT)), 1);
CPD.nsamples = 0;    
