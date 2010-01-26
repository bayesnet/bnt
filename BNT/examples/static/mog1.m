% Fit a mixture of Gaussians using netlab and BNT

rand('state', 0);
randn('state', 0);

% Q -> Y
ncenters = 2; dim = 2;
cov_type = 'full';

% Generate the data from a mixture of 2 Gaussians
%mu = randn(dim, ncenters);
mu = zeros(dim, ncenters);
mu(:,1) = [-1 -1]';
mu(:,1) = [1 1]';
Sigma = repmat(0.1*eye(dim),[1 1 ncenters]);
ndat1 = 8; ndat2 = 8;
%ndat1 = 2; ndat2 = 2;
ndata = ndat1+ndat2;
x1 = gsamp(mu(:,1), Sigma(:,:,1), ndat1);
x2 = gsamp(mu(:,2), Sigma(:,:,2), ndat2);
data = [x1; x2];
%plot(x1(:,1),x1(:,2),'ro', x2(:,1),x2(:,2),'bx')

% Fit using netlab
max_iter = 3;
mix = gmm(dim, ncenters, cov_type);
options = foptions;
options(1) = 1; % verbose
options(14) = max_iter;

% extract initial params
%mix = gmminit(mix, x, options); % Initialize with K-means
mu0 = mix.centres';
pi0 = mix.priors(:);
Sigma0 = mix.covars; % repmat(eye(dim), [1 1 ncenters]);

[mix, options] = gmmem(mix, data, options);

% Final params
ll1 = options(8);
mu1 = mix.centres';
pi1 = mix.priors(:);
Sigma1 = mix.covars;




% BNT

dag = zeros(2);
dag(1,2) = 1;
node_sizes = [ncenters dim];
discrete_nodes = 1;
onodes = 2;

bnet = mk_bnet(dag, node_sizes, 'discrete', discrete_nodes, 'observed', onodes);
bnet.CPD{1} = tabular_CPD(bnet, 1, pi0);
bnet.CPD{2} = gaussian_CPD(bnet, 2, 'mean', mu0, 'cov', Sigma0, 'cov_type', cov_type, ...
			   'cov_prior_weight', 0);

engine = jtree_inf_engine(bnet);

evidence = cell(2, ndata);
evidence(2,:) = num2cell(data', 1);

[bnet2, LL] = learn_params_em(engine, evidence, max_iter);

ll2 = LL(end);
s1 = struct(bnet2.CPD{1});
pi2 = s1.CPT(:);

s2 = struct(bnet2.CPD{2});
mu2 = s2.mean;
Sigma2 = s2.cov;

% assert(approxeq(ll1, ll2)); % gmmem returns the value after the final M step, GMT before
assert(approxeq(mu1, mu2));
assert(approxeq(Sigma1, Sigma2))
assert(approxeq(pi1, pi2))


