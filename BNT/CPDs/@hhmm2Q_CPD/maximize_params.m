function CPD = maximize_params(CPD, temp)
% MAXIMIZE_PARAMS Set the params of a hhmmQ node to their ML/MAP values.
% CPD = maximize_params(CPD, temperature)

if sum(CPD.start_counts(:)) > 0
  CPD.startprob = mk_stochastic(CPD.start_counts);
end
if sum(CPD.trans_counts(:)) > 0
  CPD.transprob = mk_stochastic(CPD.trans_counts);
end
