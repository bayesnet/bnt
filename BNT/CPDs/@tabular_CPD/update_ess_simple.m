function CPD = update_ess_simple(CPD, counts)
% UPDATE_ESS_SIMPLE Update the Expected Sufficient Statistics of a tabular node.
% function CPD = update_ess_simple(CPD, counts)

CPD.nsamples = CPD.nsamples + 1;            
CPD.counts = CPD.counts + counts(:);
