%D -> Discrete, C -> Continuous
%The model is very simple - D->C.
%This tests the ML and EM for a fully observed model using the Von Mises
%distribution as the cts distribution. The learned model parameters are
%compared to those learned from  the scipy.stats. The results can be found
%at https://github.com/sachag678/matlab_projects/blob/master/Testing%20Von%20Mises%20fit%20in%20python.ipynb

D=1;C=2;
n=2;
dag = zeros(n,n);
dag(D,C)=1;

ns = [2 1];
dnodes = 1;

bnet = mk_bnet(dag, ns, dnodes); % create bnet with dag and definition of dnodes

bnet.CPD{2} = vonMises_CPD(bnet, 2); %vonMises
bnet.CPD{1} = tabular_CPD(bnet, 1); %tabular

%get data
data = [1 0.1*pi; 2 0.8*pi; 2 0.9*pi; 1 0.2*pi];
ncases = size(data, 1);	% number of data points
cases = cell(n,ncases);	% create an empty table to store the data to be given to the learning algorithm
cases([1:n],:) = num2cell(data(:,:)');	% copy the data

%learn
%learn EM
engine = jtree_inf_engine(bnet);
[bnet2, ~, engine] = learn_params_em(engine, cases);

s = struct(bnet2.CPD{2});
mu = s.mean;
k = s.con;

%Assertion for testing the distribution parameters - against the 
%distribution params from the stats toolbox in python.  
assert(abs(mu(1)-0.4712)<0.0001)
assert(abs(mu(2)-2.6704)<0.0001)
assert(abs(k(:,:,1)-40.8667)<0.01)
assert(abs(k(:,:,2)-40.8667)<0.01)

%learn ML
bnet3 = learn_params(bnet, cases);
engine = jtree_inf_engine(bnet3);

s = struct(bnet3.CPD{2});
mu = s.mean;
k = s.con;

%Assertion for testing the distribution parameters - against the 
%distribution params from the stats toolbox in python.  
assert(abs(mu(1)-0.4712)<0.0001, 'mu(1) incorrectly calculated')
assert(abs(mu(2)-2.6704)<0.0001, 'mu(2) incorrectly calculated')
assert(abs(k(:,:,1)-40.8667)<0.01 , 'k(1) incorrectly calculated')
assert(abs(k(:,:,2)-40.8667)<0.01, 'k(1) incorrectly calculated')

%PERFORM INFERENCE
%-------------------------------------------------------------------------
evidence = cell(1,n);

evidence{C} =0.5005*pi;

engine = enter_evidence(engine, evidence);

%calculate marginal on a specific node
marg = marginal_nodes(engine, 1);
prob = marg.T;

%Check probability - calculated by hand using bayes rule.
%P(X=1,Y=0.9pi) = P(Y=0.9pi,X=1)P(X=1)/P(Y=0.9pi)
%.. = P(Y=0.9pi,X=1)P(X=1)/P(Y=0.9pi,X=1)P(X=1)+ P(Y=0.9pi,X=0)P(X=0)
%.. = 4.8877e-10/(4.8877e-10+5.4786e-10) = 0.4715
assert(abs(prob(1)-0.4715)<0.001, 'prob(1) is incorrect')
assert(abs(prob(2)-0.5285)<0.001, 'prob(1) is incorrect')

