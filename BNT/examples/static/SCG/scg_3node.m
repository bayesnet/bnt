% This example is from Page.143 of "Probabilistic Networks and Expert Systems",
% Cowell, Dawid, Lauritzen and Spiegelhalter, 1999, Springer.

X = 1; Y = 2; Z = 3;
n = 3;

dag = zeros(n);
dag(X, Y)=1;
dag(Y, Z)=1;

ns = ones(1, n);
dnodes = [];

bnet = mk_bnet(dag, ns, dnodes);
bnet.CPD{X} = gaussian_CPD(bnet, X, 'mean', 0, 'cov', 1);
bnet.CPD{Y} = gaussian_CPD(bnet, Y, 'mean', 0, 'cov', 1, 'weights', 1);
bnet.CPD{Z} = gaussian_CPD(bnet, Z, 'mean', 0, 'cov', 1, 'weights', 1);

engines = {};
engines{end+1} = jtree_inf_engine(bnet);
engines{end+1} = stab_cond_gauss_inf_engine(bnet);
nengines = length(engines);

evidence = cell(1,n);
evidence{Y} = 1.5; 

for e=1:nengines
  engines{e} = enter_evidence(engines{e}, evidence);
  margX = marginal_nodes(engines{e}, X);
  assert(approxeq(margX.mu, 0.75))
  assert(approxeq(margX.Sigma, 0.5))
  
  margZ = marginal_nodes(engines{e}, Z);
  assert(approxeq(margZ.mu, 1.5))
  assert(approxeq(margZ.Sigma, 1))
end


evidence = cell(1,n);
evidence{Z} = 1.5; 

for e=1:nengines
  engines{e} = enter_evidence(engines{e}, evidence);
  margX = marginal_nodes(engines{e}, X);
  assert(approxeq(margX.mu, 1/2))
  assert(approxeq(margX.Sigma, 2/3))
  
  margY = marginal_nodes(engines{e}, Y);
  assert(approxeq(margY.mu, 1))
  assert(approxeq(margY.Sigma, 2/3))
end

