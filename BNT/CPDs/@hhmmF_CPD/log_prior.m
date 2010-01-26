function L = log_prior(CPD)
% LOG_PRIOR Return log P(theta) for a hhmm F CPD 
% L = log_prior(CPD)

L = log_prior(CPD.sub_CPD_term);
