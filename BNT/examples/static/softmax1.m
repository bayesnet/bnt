% Check that softmax works with a simple classification demo.
% Based on netlab's demglm2
% X -> Q where X is an input node, and Q is a softmax

rand('state', 0);
randn('state', 0);
  
% Check inference

input_dim = 2;
num_classes = 3;
IRLS_iter = 3;

net = glm(input_dim, num_classes, 'softmax');

dag = zeros(2);
dag(1,2) = 1;
discrete_nodes = [2];
bnet = mk_bnet(dag, [input_dim num_classes], 'discrete', discrete_nodes, 'observed', 1);
bnet.CPD{1} = root_CPD(bnet, 1);
clamped = 0;
bnet.CPD{2} = softmax_CPD(bnet, 2, net.w1, net.b1, clamped, IRLS_iter);

engine = jtree_inf_engine(bnet);

x = rand(1, input_dim);
q = glmfwd(net, x);

[engine, ll] = enter_evidence(engine, {x, []});
m = marginal_nodes(engine, 2);
assert(approxeq(m.T(:), q(:)));


% Check learning
% We use EM, but in fact there is no hidden data.
% The M step will call IRLS on the softmax node.

% Generate data from three classes in 2d
input_dim = 2;
num_classes = 3;

% Fix seeds for reproducible results
randn('state', 42);
rand('state', 42);

ndata = 10;
% Generate mixture of three Gaussians in two dimensional space
data = randn(ndata, input_dim);
targets = zeros(ndata, 3);

% Priors for the clusters
prior(1) = 0.4;
prior(2) = 0.3;
prior(3) = 0.3;

% Cluster centres
c = [2.0, 2.0; 0.0, 0.0; 1, -1];

ndata1 = prior(1)*ndata;
ndata2 = (prior(1) + prior(2))*ndata;
% Put first cluster at (2, 2)
data(1:ndata1, 1) = data(1:ndata1, 1) * 0.5 + c(1,1);
data(1:ndata1, 2) = data(1:ndata1, 2) * 0.5 + c(1,2);
targets(1:ndata1, 1) = 1;

% Leave second cluster at (0,0)
data((ndata1 + 1):ndata2, :) = ...
  data((ndata1 + 1):ndata2, :);
targets((ndata1+1):ndata2, 2) = 1;

data((ndata2+1):ndata, 1) = data((ndata2+1):ndata,1) *0.6 + c(3, 1);
data((ndata2+1):ndata, 2) = data((ndata2+1):ndata,2) *0.6 + c(3, 2);
targets((ndata2+1):ndata, 3) = 1;


if 0
  ndata = 1;
  data = x;
  targets = [1 0 0];
end

options = foptions;
options(1) = -1; % verbose
options(14) = IRLS_iter;
[net2, options2] = glmtrain(net, options, data, targets);
net2.ll = options2(8); % type 'help foptions' for details

cases = cell(2, ndata);
for l=1:ndata
  q = find(targets(l,:)==1);
  x = data(l,:);
  cases{1,l} = x(:);
  cases{2,l} = q;
end

max_iter = 2; % we have complete observability, so 1 iter is enough
[bnet2, ll2] = learn_params_em(engine, cases, max_iter);

w = get_field(bnet2.CPD{2},'weights');
b = get_field(bnet2.CPD{2},'offset')';

w
net2.w1

b
net2.b1

% assert(approxeq(net2.ll, ll2)); % glmtrain returns ll after final M step, learn_params before

