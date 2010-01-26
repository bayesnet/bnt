% Sigmoid Belief Hierarchical Mixtures of Experts

clear all
clc
X = 1;
Q1 = 2;
Q2 = 3;
Y = 4;
dag = zeros(4,4);
dag(X,[Q1 Q2 Y]) = 1;
dag(Q1, [Q2 Y]) = 1;
dag(Q2,Y)=1;
ns = [1 3 4 3];
dnodes = [2 3 4];
onodes=[1 2 3 4];
bnet = mk_bnet(dag,ns, dnodes);

rand('state',0); randn('state',0);

bnet.CPD{1} = root_CPD(bnet, 1);
bnet.CPD{2} = softmax_CPD(bnet, 2, 'max_iter', 3);
bnet.CPD{3} = softmax_CPD(bnet, 3, 'discrete', [2], 'max_iter', 3);
bnet.CPD{4} = softmax_CPD(bnet, 4, 'discrete', [2 3], 'max_iter', 3);

T=5;
cases = cell(4, T);
cases(1,:)=num2cell(rand(1,T));
%cases(2,:)=num2cell(round(rand(1,T)*2)+1);
%cases(3,:)=num2cell(round(rand(1,T)*3)+1);
cases(4,:)=num2cell(round(rand(1,T)*2)+1);

engine = jtree_inf_engine(bnet, onodes);

[engine, loglik] = enter_evidence(engine, cases);

disp('learning-------------------------------------------')
[bnet2, LL2] = learn_params_em(engine, cases, 4);