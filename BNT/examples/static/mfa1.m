% Factor analysis
% Z -> X,  Z in R^k, X in R^D, k << D (high dimensional observations explained by small source)
% Z ~ N(0,I),   X|Z ~ N(L z, Psi), where Psi is diagonal.
%
% Mixtures of FA
% Now X|Z,W=i ~ N(mu(i) + L(i) Z, Psi(i))
%
% We compare to Zoubin Ghahramani's code.

randn('state', 0);
max_iter = 3;
M = 2;
k = 3;
D = 5;

n = 5;
X1 = randn(n, D);
X2 = randn(n, D) + 2; % move the mean to (2,2,2...)
X = [X1; X2];
N = size(X, 1);

% initialise as in mfa
tiny=exp(-700);
mX = mean(X);
cX=cov(X);
scale=det(cX)^(1/D);
randn('state',0); % must reset seed here so initial params are identical to mfa
L0=randn(D*M,k)*sqrt(scale/k);
W0 = permute(reshape(L0, [D M k]), [1 3 2]); % use D,K,M 
Psi0=diag(cX)+tiny;
Pi0=ones(M,1)/M;
Mu0=randn(M,D)*sqrtm(cX)+ones(M,1)*mX;

[Lh1, Ph1, Mu1, Pi1, LL1] = mfa(X,M,k,max_iter);
Lh1 = permute(reshape(Lh1, [D M k]), [1 3 2]); % use D,K,M 


ns = [M k D];
dag = zeros(3);
dag(1,3) = 1;
dag(2,3) = 1;
dnodes = 1;
onodes = 3;

bnet = mk_bnet(dag, ns, 'discrete', dnodes, 'observed', onodes);
bnet.CPD{1} = tabular_CPD(bnet, 1, Pi0);

%bnet.CPD{2} = gaussian_CPD(bnet, 2, zeros(k, 1), eye(k), [], 'diag', 'untied', 'clamp_mean',  'clamp_cov');

bnet.CPD{2} = gaussian_CPD(bnet, 2, 'mean', zeros(k, 1), 'cov', eye(k), 'cov_type', 'diag', ...
			   'cov_prior_weight', 0, 'clamp_mean', 1, 'clamp_cov', 1);

%bnet.CPD{3} = gaussian_CPD(bnet, 3, Mu0', repmat(diag(Psi0), [1 1 M]), W0, 'diag', 'tied');

bnet.CPD{3} = gaussian_CPD(bnet, 3, 'mean', Mu0', 'cov', repmat(diag(Psi0), [1 1 M]), ...
			   'weights', W0, 'cov_type', 'diag', 'cov_prior_weight', 0, 'tied_cov', 1);

engine = jtree_inf_engine(bnet);
evidence = cell(3, N);
evidence(3,:) = num2cell(X', 1);

[bnet2, LL2, engine2] = learn_params_em(engine, evidence, max_iter);

s = struct(bnet2.CPD{1});
Pi2 = s.CPT(:);
s = struct(bnet2.CPD{3});
Mu2 = s.mean;
W2 = s.weights;
Sigma2 = s.cov;


% Compare to Zoubin's code
assert(approxeq(LL1,LL2));
for i=1:M
  assert(approxeq(W2(:,:,i), Lh1(:,:,i)));
  assert(approxeq(Sigma2(:,:,i), diag(Ph1)));
  assert(approxeq(Mu2(:,i), Mu1(i,:)));
  assert(approxeq(Pi2(:), Pi1(:)));
end

