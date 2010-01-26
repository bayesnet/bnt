function L = log_prior(CPD)
% LOG_PRIOR Return log P(theta) for a hhmm CPD 
% L = log_prior(CPD)

L = log_prior(CPD.sub_CPD_trans);
if ~isempty(CPD.sub_CPD_start)
  L = L + log_prior(CPD.sub_CPD_start);
end
