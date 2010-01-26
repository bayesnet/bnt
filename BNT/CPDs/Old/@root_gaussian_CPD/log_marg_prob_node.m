function L = log_marg_prob_node(CPD, self_ev, pev)
% LOG_MARG_PROB_NODE Compute prod_m log P(x(i,m)| x(pi_i,m)) for node i (root_gaussian)
% L = log_marg_prob_node(CPD, self_ev, pev)
%
% This differs from log_prob_node because we integrate out the parameters.
% self_ev{m} is the evidence on this node in case m.
% pev{i,m} is the evidence on the i'th parent in case m (ignored).

ncases = length(self_ev);

if ncases==0
  L = 0;
  return;
elseif ncases==1 
  x = cat(1, self_ev{:});
  k = length(x);
  n0 = CPD.prior.n;
  mu = CPD.prior.mu;
  alpha = CPD.prior.alpha;
  beta = CPD.prior.beta;
  gamma = 2*alpha - k + 1;
  % Bernardo and Smith p441
  L = log_student_pdf(x, mu, n0/(n0+1)*0.5*gamma*inv(beta), gamma);
else
  error('can''t handle batch data');
end
