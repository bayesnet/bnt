% Check sampling on a mixture of experts model
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
ns = [1 2 2];
dnodes = [2];
bnet = mk_bnet(dag, ns, dnodes);

x = 0.5;
bnet.CPD{1} = root_CPD(bnet, 1, x);
bnet.CPD{2} = softmax_CPD(bnet, 2);
bnet.CPD{3} = gaussian_CPD(bnet, 3);

data_case = sample_bnet(bnet, 'evidence', {0.8, [], []})
ll = log_lik_complete(bnet, data_case)

data_case = sample_bnet(bnet, 'evidence', {-11, [], []})
ll = log_lik_complete(bnet, data_case)


