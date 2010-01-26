% Factor analysis
% Z -> X,  Z in R^k, X in R^D, k << D (high dimensional observations explained by small source)
% Z ~ N(0,I),   X|Z ~ N(L z, Psi), where Psi is diagonal.
%
% We compare to Zoubin Ghahramani's code.

state = 0;
rand('seed', state);
randn('seed', state);
max_iter = 3;
k = 2;
D = 4;
N = 10;
X = randn(N, D);

% Initialize as in Zoubin's ffa (fast factor analysis)
X=X-ones(N,1)*mean(X);
XX=X'*X/N;
diagXX=diag(XX);
cX=cov(X);
scale=det(cX)^(1/D);
randn('seed', 0);  % must reset seed here so initial params are identical to mfa
L0=randn(D,k)*sqrt(scale/k);
W0 = L0;
Psi0=diag(cX);

[L1, Psi1, LL1] = ffa(X,k,max_iter);


ns = [k D];
dag = zeros(2,2);
dag(1,2) = 1;
bnet = mk_bnet(dag, ns, 'discrete', [], 'observed', 2);
bnet.CPD{1} = gaussian_CPD(bnet, 1, 'mean', zeros(k,1), 'cov', eye(k), 'cov_type', 'diag', ...
			   'clamp_mean', 1, 'clamp_cov', 1);
bnet.CPD{2} = gaussian_CPD(bnet, 2, 'mean', zeros(D,1), 'cov', diag(Psi0), 'weights', W0, ...
			   'cov_type', 'diag', 'cov_prior_weight', 0, 'clamp_mean', 1);

engine = jtree_inf_engine(bnet);
evidence = cell(2,N);
evidence(2,:) = num2cell(X', 1);

[bnet2, LL2] = learn_params_em(engine, evidence, max_iter);

s = struct(bnet2.CPD{2});
L2 = s.weights;
Psi2 = s.cov;



% Compare to Zoubin's code
assert(approxeq(LL2, LL1));
assert(approxeq(Psi2, diag(Psi1)));
assert(approxeq(L2, L1));



