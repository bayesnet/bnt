function L = log_marg_prob_node(CPD, self_ev, pev)
% LOG_MARG_PROB_NODE Compute prod_m log P(x(i,m)| x(pi_i,m)) for node i (linear_gaussian)
% L = log_marg_prob_node(CPD, self_ev, pev)
%
% This differs from log_prob_node because we integrate out the parameters.
% self_ev{m} is the evidence on this node in case m.
% pev{i,m} is the evidence on the i'th parent in case m 
% We assume there is <= 1 case.

ncases = length(self_ev);

if ncases==0
  L = 0;
  return;
elseif ncases==1 
  y = self_ev{1};
  x = cat(1, pev{:}); % column vector
  f = 1-x'*inv(x*x' + CPD.prior.n)*x;
  alpha = CPD.prior.alpha;
  L = log_student_pdf(y, x'*CPD.prior.theta, f*alpha/CPD.prior.beta, 2*alpha);
else
  error('can''t handle batch data');
end
