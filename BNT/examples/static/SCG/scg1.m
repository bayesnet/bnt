% Same as cg1, except we call stab_cond_gauss_inf_engine

bnet  = mk_incinerator_bnet;

engines = {};
engines{end+1} = stab_cond_gauss_inf_engine(bnet);
engines{end+1} = jtree_inf_engine(bnet);
engines{end+1} = cond_gauss_inf_engine(bnet);
nengines = length(engines);

F = 1; W = 2; E = 3; B = 4; C = 5; D = 6; Min = 7; Mout = 8; L = 9;
n = 9;
dnodes = [B F W];
cnodes = mysetdiff(1:n, dnodes);

evidence = cell(1,n); % no evidence
ll = zeros(1, nengines);
for e=1:nengines
  [engines{e}, ll(e)] = enter_evidence(engines{e}, evidence);
end
%assert(approxeq(ll(1), ll)))
ll

% Compare to the results in table on p1107.
% These results are printed to 3dp in Cowell p150

mu = zeros(1,n);
sigma = zeros(1,n);
dprob = zeros(1,n);
addev = 1;
tol = 1e-2;
for e=1:nengines
  for i=cnodes(:)'
    m = marginal_nodes(engines{e}, i, addev);
    mu(i) = m.mu;
    sigma(i) = sqrt(m.Sigma);
  end
  for i=dnodes(:)'
    m = marginal_nodes(engines{e}, i, addev);
    dprob(i) = m.T(1);
  end
  assert(approxeq(mu([E D C L Min Mout]), [-3.25 3.04 -1.85 1.48 -0.214 2.83], tol))
  assert(approxeq(sigma([E D C L Min Mout]), [0.709 0.770 0.507 0.631 0.459 0.860], tol))
  assert(approxeq(dprob([B F W]), [0.85 0.95 0.29], tol))
  %m = marginal_nodes(engines{e}, bnet.names('E'), addev);
  %assert(approxeq(m.mu, -3.25, tol))
  %assert(approxeq(sqrt(m.Sigma), 0.709, tol))
end

% Add evidence (p 1105, top right)
evidence = cell(1,n);
evidence{W} = 1; % industrial
evidence{L} = 1.1;
evidence{C} = -0.9;

ll = zeros(1, nengines);
for e=1:nengines
  [engines{e}, ll(e)] = enter_evidence(engines{e}, evidence);
end
%assert(all(approxeq(ll(1), ll)))
ll

for e=1:nengines
  for i=cnodes(:)'
    m = marginal_nodes(engines{e}, i, addev);
    mu(i) = m.mu;
    sigma(i) = sqrt(m.Sigma);
  end
  for i=dnodes(:)'
    m = marginal_nodes(engines{e}, i, addev);
    dprob(i) = m.T(1);
  end
  assert(approxeq(mu([E D C L Min Mout]), [-3.90 3.61 -0.9 1.1 0.5 4.11], tol))
  assert(approxeq(sigma([E D C L Min Mout]), [0.076 0.326 0 0 0.1 0.344], tol))
  assert(approxeq(dprob([B F W]), [0.0122 0.9995 1], tol))
end

