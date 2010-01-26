function lam_msg = CPD_to_lambda_msg(CPD, msg_type, n, ps, msg, p, evidence)
% CPD_TO_LAMBDA_MSG Compute lambda message (gmux)
% lam_msg = compute_lambda_msg(CPD, msg_type, n, ps, msg, p, evidence)
% Pearl p183 eq 4.52

% Let Y be this node, X1..Xn be the cts parents and M the discrete switch node.
% e.g., for n=3, M=1
%
%  X1 X2 X3 M
%   \
%    \
%      Y
%
% So the only case in which we send an informative message is if p=1=M.
% To the other cts parents, we send the "know nothing" message.

switch msg_type
 case 'd',
  error('gaussian_CPD can''t create discrete msgs')
 case 'g',
  cps = ps(CPD.cps);
  cpsizes = CPD.sizes(CPD.cps);
  self_size = CPD.sizes(end);
  i = find_equiv_posns(p, cps); % p is n's i'th cts parent
  psz = cpsizes(i);
  dps = ps(CPD.dps);
  M = evidence{dps};
  if isempty(M)
    error('gmux node must have observed discrete parent')
  end
  P = msg{n}.lambda.precision;
  if all(P == 0) | (cps(M) ~= p) % if we know nothing, or are sending to a disconnected parent
    lam_msg.precision = zeros(psz, psz);
    lam_msg.info_state = zeros(psz, 1);
    return;
  end
  % We are sending a message to the only effectively connected parent.
  % There are no other incoming pi messages.
  Bmu = CPD.mean(:,M);
  BSigma = CPD.cov(:,:,M);
  Bi = CPD.weights(:,:,M);
  if (det(P) > 0) | isinf(P) 
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
    % method that uses matrix inversion lemma
    A = inv(P + inv(BSigma));
    C = P - P*A*P;
    lam_msg.precision = Bi' * C * Bi;
    D = eye(self_size) - P*A;
    z = msg{n}.lambda.info_state;
    lam_msg.info_state = Bi' * (D*z - D*P*Bmu);
  end
end
