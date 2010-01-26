function CPD = reset_ess(CPD)
% RESET_ESS Reset the Expected Sufficient Statistics for a CPD (dsoftmax)
% CPD = reset_ess(CPD)

CPD.parent_vals = [];
CPD.eso_weights=[];
CPD.self_vals = [];
CPD.nsamples = 0;  
