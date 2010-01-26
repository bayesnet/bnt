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

IRLS_iter = 10;
clamped = 0;

bnet.CPD{1} = root_CPD(bnet, 1);

% start with good initial params
w = [-5 5];  % w(:,i) is the normal vector to the i'th decisions boundary
b = [0 0];  % b(i) is the offset (bias) to the i'th decisions boundary

mu = [0 0];
sigma = 1;
Sigma = repmat(sigma*eye(ns(Y)), [ns(Y) ns(Y) ns(Q)]);
W = [-1 1];
W2 = reshape(W, [ns(Y) ns(X) ns(Q)]);

bnet.CPD{2} = softmax_CPD(bnet, 2, w, b,  clamped, IRLS_iter);
bnet.CPD{3} = gaussian_CPD(bnet, 3, 'mean', mu, 'cov', Sigma, 'weights', W2);


engine = jtree_inf_engine(bnet);

evidence = cell(1,3);
evidence{X} = 0.68;

engine = enter_evidence(engine, evidence);

m = marginal_nodes(engine, Y);
m.mu
