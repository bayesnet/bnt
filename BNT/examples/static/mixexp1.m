% Fit a piece-wise linear regression model.
% Here is the model
%
%  X \
%  | |
%  Q |
%  | /
%  Y
%
% where all arcs point down.
% We condition everything on X, so X is a root node. Q is a softmax, and Y is a linear Gaussian.
% Q is hidden, X and Y are observed.

X = 1;
Q = 2;
Y = 3;
dag = zeros(3,3);
dag(X,[Q Y]) = 1;
dag(Q,Y) = 1;
ns = [1 2 1]; % make X and Y scalars, and have 2 experts
dnodes = [2];
onodes = [1 3];
bnet = mk_bnet(dag, ns, 'discrete', dnodes, 'observed', onodes);


w = [-5 5];  % w(:,i) is the normal vector to the i'th decisions boundary
b = [0 0];  % b(i) is the offset (bias) to the i'th decisions boundary

mu = [0 0];
sigma = 1;
Sigma = repmat(sigma*eye(ns(Y)), [ns(Y) ns(Y) ns(Q)]);
W = [-1 1];
W2 = reshape(W, [ns(Y) ns(X) ns(Q)]);

bnet.CPD{1} = root_CPD(bnet, 1);
bnet.CPD{2} = softmax_CPD(bnet, 2, w, b);
bnet.CPD{3} = gaussian_CPD(bnet, 3, 'mean', mu, 'cov', Sigma, 'weights', W2);



% Check inference

x = 0.1;
ystar = 1;

engine = jtree_inf_engine(bnet);
[engine, loglik] = enter_evidence(engine, {x, [], ystar});
Qpost = marginal_nodes(engine, 2);

% eta(i,:) = softmax (gating) params for expert i
eta = [b' w'];

% theta(i,:) = regression vector for expert i
theta = [mu' W'];

% yhat(i) = E[y | Q=i, x] = prediction of i'th expert
x1 = [1 x]';
yhat = theta * x1;

% gate_prior(i,:) = Pr(Q=i | x)
gate_prior = normalise(exp(eta * x1));

% cond_lik(i) = Pr(y | Q=i, x)
cond_lik = (1/(sqrt(2*pi)*sigma)) * exp(-(0.5/sigma^2) * ((ystar - yhat) .* (ystar - yhat)));

% gate_posterior(i,:) = Pr(Q=i | x, y)
[gate_posterior, lik] = normalise(gate_prior .* cond_lik);

assert(approxeq(gate_posterior(:), Qpost.T(:)));
assert(approxeq(log(lik), loglik));


