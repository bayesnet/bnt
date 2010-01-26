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

if 0
  % start with good initial params
  w = [-5 5];  % w(:,i) is the normal vector to the i'th decisions boundary
  b = [0 0];  % b(i) is the offset (bias) to the i'th decisions boundary
  
  mu = [0 0];
  sigma = 1;
  Sigma = repmat(sigma*eye(ns(Y)), [ns(Y) ns(Y) ns(Q)]);
  W = [-1 1];
  W2 = reshape(W, [ns(Y) ns(X) ns(Q)]);

  bnet.CPD{2} = softmax_CPD(bnet, 2, w, b,  clamped, IRLS_iter);
  bnet.CPD{3} = gaussian_CPD(bnet, 3, mu, Sigma, W2);
else
  % start with rnd initial params
  rand('state', 0);
  randn('state', 0);
  bnet.CPD{2} = softmax_CPD(bnet, 2, 'clamped', clamped, 'max_iter', IRLS_iter);
  bnet.CPD{3} = gaussian_CPD(bnet, 3);
end



load('/examples/static/Misc/mixexp_data.txt', '-ascii');        
% Just use 1/10th of the data, to speed things up
data = mixexp_data(1:10:end, :);
%data = mixexp_data;
 
%plot(data(:,1), data(:,2), '.')


s = struct(bnet.CPD{2}); % violate object privacy
%eta0 = [s.glim.b1; s.glim.w1]';
eta0 = [s.glim{1}.b1; s.glim{1}.w1]';
s = struct(bnet.CPD{3}); % violate object privacy
W = reshape(s.weights, [1 2]);
theta0 = [s.mean; W]';

%figure(1)
%mixexp_plot(theta0, eta0, data);
%suptitle('before learning')

ncases = size(data, 1);
cases = cell(3, ncases);
cases([1 3], :) = num2cell(data');

engine = jtree_inf_engine(bnet);

% log lik before learning
ll = 0;
for l=1:ncases
  ev = cases(:,l);
  [engine, loglik] = enter_evidence(engine, ev);
  ll = ll + loglik;
end

% do learning
max_iter = 5;
[bnet2, LL2] = learn_params_em(engine, cases, max_iter);

s = struct(bnet2.CPD{2});
%eta2 = [s.glim.b1; s.glim.w1]';
eta2 = [s.glim{1}.b1; s.glim{1}.w1]';
s = struct(bnet2.CPD{3});
W = reshape(s.weights, [1 2]);
theta2 = [s.mean; W]';

%figure(2)
%mixexp_plot(theta2, eta2, data);
%suptitle('after learning')

fprintf('mixexp2: loglik before learning %f, after %d iters %f\n', ll, length(LL2),  LL2(end));



