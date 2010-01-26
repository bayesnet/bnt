% Sigmoid Belief IOHMM
% Here is the model
%
%  X \  X \
%  | |  | |
%  Q-|->Q-|-> ...
%  | /  | /
%  Y    Y
%
clear all;
clc;
rand('state',0); randn('state',0);
X = 1; Q = 2; Y = 3;
% intra time-slice graph
intra=zeros(3);
intra(X,[Q Y])=1;
intra(Q,Y)=1;
% inter time-slice graph
inter=zeros(3);
inter(Q,Q)=1;

ns = [1 3 1]; 
dnodes = [2];
eclass1 = [1 2 3];
eclass2 = [1 4 3];
bnet = mk_dbn(intra, inter, ns, dnodes, eclass1, eclass2);
bnet.CPD{1} = root_CPD(bnet, 1);
% ==========================================================
bnet.CPD{2} = softmax_CPD(bnet, 2);
bnet.CPD{4} = softmax_CPD(bnet, 5, 'discrete', [2]);
% ==========================================================
bnet.CPD{3} = gaussian_CPD(bnet, 3);

% make some data
T=20;
cases = cell(3, T);
cases(1,:)=num2cell(round(rand(1,T)*2)+1);
%cases(2,:)=num2cell(round(rand(1,T))+1);
cases(3,:)=num2cell(rand(1,T));

engine = bk_inf_engine(bnet, 'exact', [1 2 3]);

% log lik before learning
[engine, loglik] = enter_evidence(engine, cases);

% do learning
ev=cell(1,1);
ev{1}=cases;
[bnet2, LL2] = learn_params_dbn_em(engine, ev, 3);