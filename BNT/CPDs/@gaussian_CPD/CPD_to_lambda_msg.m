function lam_msg = CPD_to_lambda_msg(CPD, msg_type, n, ps, msg, p, evidence)
% CPD_TO_LAMBDA_MSG Compute lambda message (gaussian)
% lam_msg = compute_lambda_msg(CPD, msg_type, n, ps, msg, p, evidence)
% Pearl p183 eq 4.52

switch msg_type
 case 'd',
  error('gaussian_CPD can''t create discrete msgs')
 case 'g',
  cps = ps(CPD.cps);
  cpsizes = CPD.sizes(CPD.cps);
  self_size = CPD.sizes(end);
  i = find_equiv_posns(p, cps); % p is n's i'th cts parent
  psz = cpsizes(i);
  if all(msg{n}.lambda.precision == 0) % no info to send on
    lam_msg.precision = zeros(psz, psz);
    lam_msg.info_state = zeros(psz, 1);
    return;
  end
  [m, Q, W] = gaussian_CPD_params_given_dps(CPD, [ps n], evidence);
  Bmu = m;
  BSigma = Q;
  for k=1:length(cps) % only get pi msgs from cts parents
    pk = cps(k);
    if pk ~= p
      %bk = block(k, cpsizes);
      bk = CPD.cps_block_ndx{k};
      Bk = W(:, bk);
      m = msg{n}.pi_from_parent{k}; 
      BSigma = BSigma + Bk * m.Sigma * Bk';
      Bmu = Bmu + Bk * m.mu;
    end
  end
  % BSigma = Q + sum_{k \neq i} B_k Sigma_k B_k'
  %bi = block(i, cpsizes);
  bi = CPD.cps_block_ndx{i};
  Bi = W(:,bi);
  P = msg{n}.lambda.precision;
  if (rcond(P) > 1e-3) | isinf(P)
    if isinf(P) % Y is observed
      Sigma_lambda = zeros(self_size, self_size); % infinite precision => 0 variance
      mu_lambda = msg{n}.lambda.mu; % observed_value;
    else
      Sigma_lambda = inv(P);
      mu_lambda = Sigma_lambda * msg{n}.lambda.info_state;
    end
    C = inv(Sigma_lambda + BSigma);
    lam_msg.precision = Bi' * C * Bi;
    lam_msg.info_state = Bi' * C * (mu_lambda - Bmu);
  else
    % method that uses matrix inversion lemma to avoid inverting P
    A = inv(P + inv(BSigma));
    C = P - P*A*P;
    lam_msg.precision = Bi' * C * Bi;
    D = eye(self_size) - P*A;
    z = msg{n}.lambda.info_state;
    lam_msg.info_state = Bi' * (D*z - D*P*Bmu);
  end
end
