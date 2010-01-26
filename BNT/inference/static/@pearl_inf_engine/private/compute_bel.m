function bel = compute_bel(msg_type, pi, lambda)

switch msg_type,
 case 'd', bel = normalise(pi .* lambda);
 case 'g',
  if isinf(lambda.precision) % ignore pi because lambda is completely certain (observed)
    bel.mu = lambda.mu;
    bel.Sigma = zeros(length(bel.mu)); % infinite precision => 0 variance
  elseif all(pi.Sigma==0) % ignore lambda because pi is completely certain (delta fn prior)
    bel.Sigma = pi.Sigma;
    bel.mu = pi.mu;
  elseif all(isinf(pi.Sigma)) % ignore pi because pi is completely uncertain
    bel.Sigma  = inv(lambda.precision);
    bel.mu = bel.Sigma * lambda.info_state;
  elseif all(lambda.precision == 0) % ignore lambda because lambda is completely uncertain
    bel.Sigma = pi.Sigma;
    bel.mu = pi.mu;
  else % combine both pi and lambda
    pi_precision = inv(pi.Sigma);
    bel.Sigma = inv(pi_precision + lambda.precision);
    bel.mu = bel.Sigma*(pi_precision * pi.mu + lambda.info_state);
  end
 otherwise, error(['unrecognized msg type ' msg_type])
end
