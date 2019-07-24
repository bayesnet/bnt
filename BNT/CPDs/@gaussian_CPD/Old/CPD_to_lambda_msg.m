function lam_msg = CPD_to_lambda_msg(CPD, msg_type, n, ps, msg, p)
% CPD_TO_LAMBDA_MSG Compute lambda message (gaussian)
% lam_msg = compute_lambda_msg(CPD, msg_type, n, ps, msg, p)
% Pearl p183 eq 4.52

switch msg_type
 case 'd',
  error('gaussian_CPD can''t create discrete msgs')
 case 'g',
  self_size = CPD.sizes(end);
  if all(msg{n}.lambda.precision == 0) % no info to send on
    lam_msg.precision = zeros(self_size);
    lam_msg.info_state = zeros(self_size, 1);
    return;
  end
  cpsizes = CPD.sizes(CPD.cps);
  dpval = 1;
  Q = CPD.cov(:,:,dpval);
  Sigmai = Q;
  wmu = zeros(self_size, 1);
  for k=1:length(ps)
    pk = ps(k);
    if pk ~= p
      bk = block(k, cpsizes);
      Bk = CPD.weights(:, bk, dpval);
      m = msg{n}.pi_from_parent{k};
      Sigmai = Sigmai + Bk * m.Sigma * Bk';
      wmu = wmu + Bk * m.mu; % m.mu = u(k)
    end
  end
  % Sigmai = Q + sum_{k \neq i} B_k Sigma_k B_k'
  i = find_equiv_posns(p, ps);
  bi = block(i, cpsizes);
  Bi = CPD.weights(:,bi, dpval);
  
  if 0
  P = msg{n}.lambda.precision;
  if isinf(P) % inv(P)=Sigma_lambda=0
    precision_temp = inv(Sigmai);
    lam_msg.precision = Bi' * precision_temp * Bi;
    lam_msg.info_state = precision_temp * (msg{n}.lambda.mu - wmu);
  else
    A = inv(P + inv(Sigmai));
    precision_temp = P + P*A*P;
    lam_msg.precision = Bi' * precision_temp * Bi;
    self_size = length(P);
    C = eye(self_size) + P*A;
    z = msg{n}.lambda.info_state;
    lam_msg.info_state = C*z - C*P*wmu;
  end
  end
  
  if isinf(msg{n}.lambda.precision)
    Sigma_lambda = zeros(self_size, self_size); % infinite precision => 0 variance
    mu_lambda = msg{n}.lambda.mu; % observed_value;
  else
    Sigma_lambda = inv(msg{n}.lambda.precision);
    mu_lambda = Sigma_lambda * msg{n}.lambda.info_state;
  end
  precision_temp = inv(Sigma_lambda + Sigmai);
  lam_msg.precision = Bi' * precision_temp * Bi;
  lam_msg.info_state = Bi' * precision_temp * (mu_lambda - wmu);
end

