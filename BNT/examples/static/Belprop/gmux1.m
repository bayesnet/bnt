% Test gmux.
% The following model, where Y is a gmux node,
% and M is set to 1, should be equivalent to X1 -> Y
%
% X1 Xn M
% \ |  /
%   Y

n = 3;
N = n+2;
Xs = 1:n;
M = n+1; 
Y = n+2;
dag = zeros(N,N);
dag([Xs M], Y)=1; 

dnodes = M;
ns = zeros(1, N);
sz = 2;
ns(Xs) = sz;
ns(M) = n;
ns(Y) = sz;

bnet = mk_bnet(dag, ns, 'discrete', M, 'observed', [M Y]);

psz = ns(Xs(1));
selfsz = ns(Y);

W = randn(selfsz, psz);
mu = randn(selfsz, 1);
Sigma = eye(selfsz, selfsz);

bnet.CPD{M} = root_CPD(bnet, M);
for i=Xs(:)'
  bnet.CPD{i} = gaussian_CPD(bnet, i, 'mean', zeros(psz, 1), 'cov', eye(psz, psz));
end
bnet.CPD{Y} = gmux_CPD(bnet, Y, 'mean', mu, 'weights', W, 'cov', Sigma);
  
evidence = cell(1,N);
yval = randn(selfsz, 1);
evidence{Y} = yval;
m = 2;
%notm = not(m-1)+1; % only valid for n=2
notm = mysetdiff(1:n, m);
evidence{M} = m;

engines = {};
engines{end+1} = jtree_inf_engine(bnet);
engines{end+1} = pearl_inf_engine(bnet, 'protocol', 'parallel');

for e=1:length(engines)
  engines{e} = enter_evidence(engines{e}, evidence);
  mXm{e} = marginal_nodes(engines{e}, Xs(m));

  % Since M=m, only Xm was updated.
  % Hence the posterior on Xnotm should equal the prior.
  for i=notm(:)'
    mXnotm = marginal_nodes(engines{e}, Xs(i));
    assert(approxeq(mXnotm.mu, zeros(psz,1)))
    assert(approxeq(mXnotm.Sigma, eye(psz, psz)))
  end
end

% Check that all engines give the same posterior
for e=2:length(engines)
  assert(approxeq(mXm{e}.mu, mXm{1}.mu))
  assert(approxeq(mXm{e}.Sigma, mXm{1}.Sigma))
end


% Compute the correct posterior by building Xm -> Y

N = 2;
dag = zeros(N,N);
dag(1, 2)=1;
ns = [psz selfsz];
bnet = mk_bnet(dag, ns, 'discrete', [], 'observed', 2);

bnet.CPD{1} = gaussian_CPD(bnet, 1, 'mean', zeros(psz, 1), 'cov', eye(psz, psz));
bnet.CPD{2} = gaussian_CPD(bnet, 2, 'mean', mu, 'cov', Sigma, 'weights', W);

jengine  = jtree_inf_engine(bnet);
evidence = {[], yval};
jengine = enter_evidence(jengine, evidence); % apply Bayes rule to invert the arc
mX = marginal_nodes(jengine, 1);

for e=1:length(engines)
  assert(approxeq(mX.mu, mXm{e}.mu))
  assert(approxeq(mX.Sigma, mXm{e}.Sigma))
end
