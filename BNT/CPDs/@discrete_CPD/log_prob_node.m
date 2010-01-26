function L = log_prob_node(CPD, self_ev, pev)
% LOG_PROB_NODE Compute sum_m log P(x(i,m)| x(pi_i,m), theta_i) for node i (discrete)
% L = log_prob_node(CPD, self_ev, pev)
%
% self_ev(m) is the evidence on this node in case m.
% pev(i,m) is the evidence on the i'th parent in case m (if there are any parents).
% (These may also be cell arrays.)

[P, p] = prob_node(CPD, self_ev, pev); % P may underflow, so we use p
tiny = exp(-700);
p = p + (p==0)*tiny; % replace 0s by tiny
L = sum(log(p));
