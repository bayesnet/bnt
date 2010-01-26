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
bnet = mk_dbn(intra, inter, ns, 'discrete', [], 'observed', 2);

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
engine{end+1} = smoother_engine(jtree_2TBN_inf_engine(bnet));
N = length(engine);


inf_time = cmp_inference_dbn(bnet, engine, T);

ncases = 2;
max_iter = 2;
[learning_time, CPD, LL, cases] = cmp_learning_dbn(bnet, engine, T, 'ncases', ncases, 'max_iter', max_iter);


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
assert(approxeq(LL2trace, LL{1}))

