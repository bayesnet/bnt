function CPD = reset_ess(CPD)
% RESET_ESS Reset the Expected Sufficient Statistics for a CPD (mlp)
% CPD = reset_ess(CPD)

CPD.W1 = [];
CPD.W2 = [];
CPD.b1 = [];
CPD.b2 = [];
CPD.parent_vals = [];
CPD.eso_weights=[];
CPD.self_vals = [];
CPD.nsamples = 0;  