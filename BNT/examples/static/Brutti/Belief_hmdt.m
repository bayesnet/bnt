% Sigmoid Belief Hidden Markov Decision Tree    (Jordan/Gharhamani 1996)
% 
clear all;
%clc;
rand('state',0); randn('state',0);
X = 1; Q1 = 2; Q2 = 3; Y = 4;
% intra time-slice graph
intra=zeros(4);
intra(X,[Q1 Q2 Y])=1;
intra(Q1,[Q2 Y])=1;
intra(Q2, Y)=1;
% inter time-slice graph
inter=zeros(4);
inter(Q1,Q1)=1;
inter(Q2,Q2)=1;

ns = [1 2 3 1]; 
dnodes = [2 3]; 
eclass1 = [1 2 3 4];
eclass2 = [1 5 6 4];
bnet = mk_dbn(intra, inter, ns, dnodes, eclass1, eclass2);

bnet.CPD{1} = root_CPD(bnet, 1);
% =========================================
bnet.CPD{2} = softmax_CPD(bnet, 2);
bnet.CPD{3} = softmax_CPD(bnet, 3, 'discrete', [2]);
bnet.CPD{5} = softmax_CPD(bnet, 6);
bnet.CPD{6} = softmax_CPD(bnet, 7, 'discrete', [3 6]);
% =========================================
bnet.CPD{4} = gaussian_CPD(bnet, 4);

% make some data
T=20;
cases = cell(4, T);
cases(1,:)=num2cell(round(rand(1,T)*2)+1);
%cases(2,:)=num2cell(round(rand(1,T))+1);
%cases(3,:)=num2cell(round(rand(1,T)*2)+1);
cases(4,:)=num2cell(rand(1,T));

engine = bk_inf_engine(bnet, 'exact', [1 2 3 4]);

% log lik before learning
[engine, loglik] = enter_evidence(engine, cases);

% do learning
ev=cell(1,1);
ev{1}=cases;
[bnet2, LL2] = learn_params_dbn_em(engine, ev, 10);