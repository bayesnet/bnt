function [h, hdata] = glmhess_weighted(net, x, t, eso_w, hdata)
%GLMHESS Evaluate the Hessian matrix for a generalised linear model.
%
%	Description
%	H = GLMHESS(NET, X, T) takes a GLM network data structure NET, a
%	matrix X of input values, and a matrix T of target values and returns
%	the full Hessian matrix H corresponding to the second derivatives of
%	the negative log posterior distribution, evaluated for the current
%	weight and bias values as defined by NET. Note that the target data
%	is not required in the calculation, but is included to make the
%	interface uniform with NETHESS.  For linear and logistic outputs, the
%	computation is very simple and is  done (in effect) in one line in
%	GLMTRAIN.
%
%	See also
%	GLM, GLMTRAIN, HESSCHEK, NETHESS
%
%	Copyright (c) Ian T Nabney (1996-9)

% Check arguments for consistency
errstring = consist(net, 'glm', x, t);
if ~isempty(errstring);
  error(errstring);
end

ndata = size(x, 1);
nparams = net.nwts;
nout = net.nout;
p = glmfwd(net, x);
inputs = [x ones(ndata, 1)];

if nargin == 4
   hdata = zeros(nparams);	% Full Hessian matrix
   % Calculate data component of Hessian
   switch net.outfn
    
   case 'softmax'
    bb_start = nparams - nout + 1;	% Start of bias weights block
    ex_hess = zeros(nparams);	% Contribution to Hessian from single example
    for m = 1:ndata
      X = x(m,:)'*x(m,:);
      a = diag(p(m,:))-((p(m,:)')*p(m,:)); 
      a=eso_w(m,1)*a;
      ex_hess(1:nparams-nout,1:nparams-nout) = kron(a, X);
      ex_hess(bb_start:nparams, bb_start:nparams) = a.*ones(net.nout, net.nout);
      temp = kron(a, x(m,:));
      ex_hess(bb_start:nparams, 1:nparams-nout) = temp;
      ex_hess(1:nparams-nout, bb_start:nparams) = temp';
      hdata = hdata + ex_hess;
    end
    
    otherwise
      error(['Unknown activation function ', net.actfn]);
    end
end

[h, hdata] = hbayes(net, hdata);
