echo on;

clc;

% This is a very basic demo of the mixture of factor analyzer software
% written in Matlab by	Zoubin Ghahramani
%			Dept of Computer Science
%			University of Toronto

pause;		% Hit any key to continue 

% To demonstrate the software we generate a sample data set
% from a mixture of two Gaussians

pause;		% Hit any key to continue 

X1=randn(300,5);	% zero mean 5 dim Gaussian data 
X2=randn(200,5)+2;	% 5 dim Gaussian data with mean [1 1 1 1 1]
X=[X1;X2];		% total 500 data points from mixture

% Fitting the model is very easy. For example to fit a mixture of 2
% factor analyzers with three factors each...

pause;		% Hit any key to continue 


[Lh,Ph,Mu,Pi,LL]=mfa(X,2,3);

% Lh, Ph, Mu, and Pi are the factor loadings, observervation
% variances, observation means for each mixture, and mixing
% proportions. LL is the vector of log likelihoods (the learning
% curve). For more information type: help mfa

% to plot the learning curve (log likelihood at each step of EM)...

pause;		% Hit any key to continue 

plot(LL);

% you get a more informative picture of convergence by looking at the
% log of the first difference of the log likelihoods...

pause;		% Hit any key to continue 

semilogy(diff(LL)); 

% you can look at some of the parameters of the fitted model... 

pause;		% Hit any key to continue 

Mu

Pi

% ...to see whether they make any sense given that me know how the
% data was generated. 

% you can also evaluate the log likelihood of another data set under
% the model we have just fitted using the mfa_cl (for Calculate
% Likelihood) function. For example, here we generate a test from the
% same distribution. 


X1=randn(300,5);
X2=randn(200,5)+2;
Xtest=[X1; X2];

pause;		% Hit any key to continue 

mfa_cl(Xtest,Lh,Ph,Mu,Pi)

% we should expect the log likelihood of the test set to be lower than
% that of the training set.

% finally, we can also fit a regular factor analyzer using the ffa
% function (Fast Factor Analysis)...

pause;		% Hit any key to continue 

[L,Ph,LL]=ffa(X,3);
  
