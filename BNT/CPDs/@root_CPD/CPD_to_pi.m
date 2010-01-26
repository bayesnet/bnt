function pi = CPD_to_pi(CPD, msg_type, n, ps, msg, evidence)
% CPD_TO_PI Compute the pi vector (root)
% function pi = CPD_to_pi(CPD, msg_type, n, ps, msg, evidence)

self_ev = evidence{n};
switch msg_type
 case 'd',
  error('root_CPD can''t create discrete msgs')
 case 'g',
  pi.mu = self_ev;
  pi.Sigma = zeros(size(self_ev));
end
