function pot = linear_vonMises_to_cpot(mu, k, W, domain, ns, cnodes, evidence)
% LINEAR_VONMISES_TO_CPOT Convert a linear Von Mises CPD  to a canonical potential.
% pot = linear_vonMises_to_cpot(mu, Sigma, W, domain, ns, cnodes, evidence)
%
% We include any cts evidence, but ignore any discrete evidence.
% (Use vonMises_CPD_params_given_dps to use discrete evidence to select mu, Sigma, W.)

odom = domain(~isemptycell(evidence(domain)));
hdom = domain(isemptycell(evidence(domain))); 
cobs = myintersect(cnodes, odom);
chid = myintersect(cnodes, hdom);
cvals = cat(1, evidence{cobs});

%[g,h,K] = vonMises_to_canonical(mu, Sigma, W);
g = log(1/(2*pi*besseli(0,k)));
if isempty(W) | (size(W,2)==0)  % no cts parents
  h = k*cos(mu);
  K = k*sin(mu);
else
  %currently not implemented
   error('cts parents not supported yet');
end

if ~isempty(cvals)
  %[g, h, K] = enter_evidence_canonical(g, h, K, chid, cobs, cvals(:), ns);  
  [hx, hy, KXX, KXY, KYX, KYY] = partition_matrix_vec(h, K, chid, cobs, ns);
  y = cvals(:);
  g = g + hy'*cos(y) + KYY*sin(y);
  if length(hx)==0 % isempty(X) % i.e., we have instantiated everything away
    h = [];
    K = [];
  else
   error('cts parents not supported yet');
  end
end

ns(odom) = 0;
pot = cpot(domain, ns(domain), g, h, K);


