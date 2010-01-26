function pi = CPD_to_pi(CPD, msg_type, n, ps, msg, evidence)
% CPD_TO_PI Compute the pi vector (gaussian)
% function pi = CPD_to_pi(CPD, msg_type, n, ps, msg, evidence)

switch msg_type
 case 'd',
  error('gaussian_CPD can''t create discrete msgs')
 case 'g',
  dps = ps(CPD.dps);
  k = evidence{dps};
  if isempty(k)
    error('gmux node must have observed discrete parent')
  end
  m = msg{n}.pi_from_parent{k}; 
  B = CPD.weights(:,:,k);
  pi.mu = CPD.mean(:,k) + B * m.mu;
  pi.Sigma = CPD.cov(:,:,k) + B * m.Sigma * B';
end
