function [e, edata, eprior, y, a] = glmerr_weighted(net, x, t, eso_w)
%GLMERR	Evaluate error function for generalized linear model.
%
%	Description
%	 E = GLMERR(NET, X, T) takes a generalized linear model data
%	structure NET together with a matrix X of input vectors and a matrix
%	T of target vectors, and evaluates the error function E. The choice
%	of error function corresponds to the output unit activation function.
%	Each row of X corresponds to one input vector and each row of T
%	corresponds to one target vector.
%
%	[E, EDATA, EPRIOR, Y, A] = GLMERR(NET, X, T) also returns the data
%	and prior components of the total error.
%
%	[E, EDATA, EPRIOR, Y, A] = GLMERR(NET, X) also returns a matrix Y
%	giving the outputs of the models and a matrix A  giving the summed
%	inputs to each output unit, where each row corresponds to one
%	pattern.
%
%	See also
%	GLM, GLMPAK, GLMUNPAK, GLMFWD, GLMGRAD, GLMTRAIN
%

%	Copyright (c) Ian T Nabney (1996-9)

% Check arguments for consistency
errstring = consist(net, 'glm', x, t);
if ~isempty(errstring);
  error(errstring);
end

[y, a] = glmfwd(net, x);

%switch net.actfn
  switch net.outfn

  case 'softmax'   	% Softmax outputs
  
    nout = size(a,2);
    % Ensure that sum(exp(a), 2) does not overflow
    maxcut = log(realmax) - log(nout);
    % Ensure that exp(a) > 0
    mincut = log(realmin);
    a = min(a, maxcut);
    a = max(a, mincut);
    temp = exp(a);
    y = temp./(sum(temp, 2)*ones(1,nout));
    % Ensure that log(y) is computable
    y(y<realmin) = realmin;
    e_app=sum(t.*log(y),2);
    edata = - eso_w'*e_app;
    
  otherwise
    error(['Unknown activation function ', net.actfn]);
end

[e, edata, eprior] = errbayes(net, edata);
