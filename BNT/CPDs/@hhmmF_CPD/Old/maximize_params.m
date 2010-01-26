function CPD = maximize_params(CPD, temp)
% MAXIMIZE_PARAMS Set the params of a hhmmF node to their ML/MAP values.
% CPD = maximize_params(CPD, temperature)

CPD.sub_CPD_term = maximize_params(CPD.sub_CPD_term, temp);
S = struct(CPD.sub_CPD_term);
CPD.termprob = S.CPT;

CPD = update_CPT(CPD);
