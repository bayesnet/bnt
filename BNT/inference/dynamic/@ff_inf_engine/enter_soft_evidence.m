function [marginals, fwd, back, loglik] = enter_soft_evidence(engine, CPDpot, observed, pot_type, filter)
% ENTER_SOFT_EVIDENCE Add the specified soft evidence to the network (ff)
% [marginals, loglik] = enter_soft_evidence(engine, CPDpot, observed, pot_type, filter)

if filter
  [fwd, loglik] = filter_evidence(engine, CPDpot, observed, pot_type);
  marginals = fwd;
  back = [];
else
  [marginals, fwd, back, loglik] = smooth_evidence(engine, CPDpot, observed, pot_type);
end
