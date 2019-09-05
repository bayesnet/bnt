function pi = CPD_to_pi(CPD, msg_type, n, ps, msg, evidence)
% CPD_TO_PI Compute the pi vector (gaussian)
% function pi = CPD_to_pi(CPD, msg_type, n, ps, msg, evidence)

switch msg_type
 case 'd',
  error('gaussian_CPD can''t create discrete msgs')
 case 'g',
  [m, Q, W] = gaussian_CPD_params_given_dps(CPD, [ps n], evidence);
  cps = ps(CPD.cps);
  cpsizes = CPD.sizes(CPD.cps);
  pi.mu = m;
  pi.Sigma = Q;
  for k=1:length(cps) % only get pi msgs from cts parents
    %bk = block(k, cpsizes);
    bk = CPD.cps_block_ndx{k};
    Bk = W(:, bk);
    m = msg{n}.pi_from_parent{k}; 
    pi.Sigma = pi.Sigma + Bk * m.Sigma * Bk';
    pi.mu = pi.mu + Bk * m.mu; % m.mu = u(k)
  end
end
