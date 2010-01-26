function pot = linear_gaussian_to_cpot(mu, Sigma, W, domain, ns, cnodes, evidence)
% LINEAR_GAUSSIAN_TO_CPOT Convert a linear Gaussian CPD  to a canonical potential.
% pot = linear_gaussian_to_cpot(mu, Sigma, W, domain, ns, cnodes, evidence)
%
% We include any cts evidence, but ignore any discrete evidence.
% (Use gaussian_CPD_params_given_dps to use discrete evidence to select mu, Sigma, W.)

odom = domain(~isemptycell(evidence(domain)));
hdom = domain(isemptycell(evidence(domain)));
cobs = myintersect(cnodes, odom);
chid = myintersect(cnodes, hdom);
cvals = cat(1, evidence{cobs});

%[g,h,K] = gaussian_to_canonical(mu, Sigma, W);
Sinv = inv(Sigma);
g = -0.5*mu'*Sinv*mu + log(normal_coef(Sigma));
if isempty(W) | (size(W,2)==0)  % no cts parents
  h = Sinv*mu;
  K = Sinv;
else
  h = [-W'*Sinv*mu; Sinv*mu];
  K = [W'*Sinv*W  -W'*Sinv';
       -Sinv*W     Sinv]; 
end

if ~isempty(cvals)
  %[g, h, K] = enter_evidence_canonical(g, h, K, chid, cobs, cvals(:), ns);  
  [hx, hy, KXX, KXY, KYX, KYY] = partition_matrix_vec(h, K, chid, cobs, ns);
  y = cvals(:);
  g = g + hy'*y - 0.5*y'*KYY*y;
  if length(hx)==0 % isempty(X) % i.e., we have instantiated everything away
    h = [];
    K = [];
  else
    h = hx - KXY*y;
    K = KXX;
  end
end

ns(odom) = 0;
pot = cpot(domain, ns(domain), g, h, K);


