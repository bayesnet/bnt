% Sigmoid Belief Net

clear all
clc
dum1 = 1;
dum2 = 2;
dum3 = 3;
Q1 = 4;
Q2 = 5;
Y = 6;
dag = zeros(6,6);
dag(dum1,[Q1 Y]) = 1;
dag(dum2, Q2)=1;
dag(dum3, [Q1 Q2])=1;
dag(Q1,[Q2 Y]) = 1;
dag(Q2, Y)=1;

ns = [2 2 3 3 4 3];
dnodes = [1:6];
bnet = mk_bnet(dag,ns, dnodes);

rand('state',0); randn('state',0);
n_iter=10;
clamped=0;

bnet.CPD{1} = tabular_CPD(bnet, 1);
bnet.CPD{2} = tabular_CPD(bnet, 2);
bnet.CPD{3} = tabular_CPD(bnet, 3);
% CPD = dsoftmax_CPD(bnet, self, dummy_pars, w, b, clamped, max_iter, verbose, wthresh,...
%    llthresh, approx_hess)
bnet.CPD{4} = softmax_CPD(bnet, 4, 'discrete', [1 3]);
bnet.CPD{5} = softmax_CPD(bnet, 5, 'discrete', [2 3]);
bnet.CPD{6} = softmax_CPD(bnet, 6, 'discrete', [1 4]);

T=5;
cases = cell(6, T);
cases(1,:)=num2cell(round(rand(1,T)*1)+1);
%cases(2,:)=num2cell(round(rand(1,T)*1)+1);
cases(3,:)=num2cell(round(rand(1,T)*2)+1);
cases(4,:)=num2cell(round(rand(1,T)*2)+1); 
%cases(5,:)=num2cell(round(rand(1,T)*3)+1);
cases(6,:)=num2cell(round(rand(1,T)*2)+1);

engine = jtree_inf_engine(bnet);

[engine, loglik] = enter_evidence(engine, cases);

disp('learning-------------------------------------------')
[bnet2, LL2, eng2] = learn_params_em(engine, cases, n_iter);