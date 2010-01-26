% to test whether scg inference engine can handl dynameic BN
% Make a linear dynamical system
%   X1 -> X2
%   |     | 
%   v     v
%   Y1    Y2 

intra = zeros(2);
intra(1,2) = 1;
inter = zeros(2);
inter(1,1) = 1;
n = 2;

X = 2; % size of hidden state
Y = 2; % size of observable state

ns = [X Y];
dnodes = [];
onodes = [2];
eclass1 = [1 2];
eclass2 = [3 2];
bnet = mk_dbn(intra, inter, ns, dnodes, eclass1, eclass2);

x0 = rand(X,1);
V0 = eye(X);
C0 = rand(Y,X);
R0 = eye(Y);
A0 = rand(X,X);
Q0 = eye(X);

bnet.CPD{1} = gaussian_CPD(bnet, 1, 'mean', x0, 'cov', V0);
%bnet.CPD{2} = gaussian_CPD(bnet, 2, 'mean', zeros(Y,1), 'cov', R0, 'weights', C0, 'full', 'untied', 'clamped_mean');
%bnet.CPD{3} = gaussian_CPD(bnet, 3, 'mean', zeros(X,1), 'cov', Q0, 'weights', A0, 'full', 'untied', 'clamped_mean');
bnet.CPD{2} = gaussian_CPD(bnet, 2, 'mean', zeros(Y,1), 'cov', R0, 'weights', C0);
bnet.CPD{3} = gaussian_CPD(bnet, 3, 'mean', zeros(X,1), 'cov', Q0, 'weights', A0);


T = 5; % fixed length sequences

clear engine;
%engine{1} = kalman_inf_engine(bnet, onodes);
engine{1} = scg_unrolled_dbn_inf_engine(bnet, T, onodes);
engine{2} = jtree_unrolled_dbn_inf_engine(bnet, T);

N = length(engine);

% inference

ev = sample_dbn(bnet, T);
evidence = cell(n,T);
evidence(onodes,:) = ev(onodes, :);

t = 2;
query = [1 3];
m = cell(1, N);
ll = zeros(1, N);

engine{1} = enter_evidence(engine{1}, evidence);
[engine{2}, ll(2)] = enter_evidence(engine{2}, evidence);
m{1} = marginal_nodes(engine{1}, query);
m{2} = marginal_nodes(engine{2}, query, t);


% compare all engines to engine{1}
for i=2:N
  assert(approxeq(m{1}.mu, m{i}.mu));
  assert(approxeq(m{1}.Sigma, m{i}.Sigma));
%  assert(approxeq(ll(1), ll(i)));
end

