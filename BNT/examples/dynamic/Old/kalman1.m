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
bnet = mk_dbn(intra, inter, ns, 'discrete', dnodes, 'eclass1', eclass1, 'eclass2', eclass2, ...
	      'observed', onodes);

x0 = rand(X,1);
V0 = eye(X);
C0 = rand(Y,X);
R0 = eye(Y);
A0 = rand(X,X);
Q0 = eye(X);

bnet.CPD{1} = gaussian_CPD(bnet, 1, 'mean', x0, 'cov', V0, 'cov_prior_weight', 0);
bnet.CPD{2} = gaussian_CPD(bnet, 2, 'mean', zeros(Y,1), 'cov', R0, 'weights', C0, ...
			   'clamp_mean', 1, 'cov_prior_weight', 0);
bnet.CPD{3} = gaussian_CPD(bnet, 3, 'mean', zeros(X,1), 'cov', Q0, 'weights', A0, ...
			   'clamp_mean', 1, 'cov_prior_weight', 0);


T = 5; % fixed length sequences

clear engine;
engine{1} = kalman_inf_engine(bnet);
engine{2} = jtree_unrolled_dbn_inf_engine(bnet, T);
engine{3} = jtree_dbn_inf_engine(bnet);
N = length(engine);

% inference

ev = sample_dbn(bnet, T);
evidence = cell(n,T);
evidence(onodes,:) = ev(onodes, :);

t = 1;
query = [1 3];
m = cell(1, N);
ll = zeros(1, N);
for i=1:N
  [engine{i}, ll(i)] = enter_evidence(engine{i}, evidence);
  m{i} = marginal_nodes(engine{i}, query, t);
end

% compare all engines to engine{1}
for i=2:N
  assert(approxeq(m{1}.mu, m{i}.mu));
  assert(approxeq(m{1}.Sigma, m{i}.Sigma));
  assert(approxeq(ll(1), ll(i)));
end

if 0
for i=2:N
  approxeq(m{1}.mu, m{i}.mu)
  approxeq(m{1}.Sigma, m{i}.Sigma)
  approxeq(ll(1), ll(i))
end
end

% learning

ncases = 5;
cases = cell(1, ncases);
for i=1:ncases
  ev = sample_dbn(bnet, T);
  cases{i} = cell(n,T);
  cases{i}(onodes,:) = ev(onodes, :);
end

max_iter = 2;
bnet2 = cell(1,N);
LLtrace = cell(1,N);
for i=1:N
  [bnet2{i}, LLtrace{i}] = learn_params_dbn_em(engine{i}, cases, 'max_iter', max_iter);
end

for i=1:N
  temp = bnet2{i};
  for e=1:3
    CPD{i,e} = struct(temp.CPD{e});
  end
end

for i=2:N
  assert(approxeq(LLtrace{i}, LLtrace{1}));
  for e=1:3
    assert(approxeq(CPD{i,e}.mean, CPD{1,e}.mean));
    assert(approxeq(CPD{i,e}.cov, CPD{1,e}.cov));
    assert(approxeq(CPD{i,e}.weights, CPD{1,e}.weights));
  end
end


% Compare to KF toolbox

data = zeros(Y, T, ncases);
for i=1:ncases
  data(:,:,i) = cell2num(cases{i}(onodes, :));
end   
[A2, C2, Q2, R2, x2, V2, LL2trace] =  learn_kalman(data, A0, C0, Q0, R0, x0, V0, max_iter);


e = 1;
assert(approxeq(x2, CPD{e,1}.mean))
assert(approxeq(V2, CPD{e,1}.cov))
assert(approxeq(C2, CPD{e,2}.weights))
assert(approxeq(R2, CPD{e,2}.cov));
assert(approxeq(A2, CPD{e,3}.weights))
assert(approxeq(Q2, CPD{e,3}.cov));
assert(approxeq(LL2trace, LLtrace{1}))

