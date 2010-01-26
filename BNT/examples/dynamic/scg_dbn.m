% Test whether stable conditional Gaussian inference works
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

bnet.CPD{1} = gaussian_CPD(bnet, 1, 'mean', x0, 'cov', V0);
bnet.CPD{2} = gaussian_CPD(bnet, 2, 'mean', zeros(Y,1), 'cov', R0, 'weights', C0);
bnet.CPD{3} = gaussian_CPD(bnet, 3, 'mean', zeros(X,1), 'cov', Q0, 'weights', A0);


T = 5; % fixed length sequences

engine = {};
engine{end+1} = kalman_inf_engine(bnet);
engine{end+1} = scg_unrolled_dbn_inf_engine(bnet, T);
engine{end+1} = jtree_unrolled_dbn_inf_engine(bnet, T);

inf_time = cmp_inference_dbn(bnet, engine, T, 'check_ll', 0);
