function [net, options] = glmtrain_weighted(net, options, x, t, eso_w, alfa)
%GLMTRAIN Specialised training of generalized linear model
%
%	Description
%	NET = GLMTRAIN(NET, OPTIONS, X, T) uses the iterative reweighted
%	least squares (IRLS) algorithm to set the weights in the generalized
%	linear model structure NET.  This is a more efficient alternative to
%	using GLMERR and GLMGRAD and a non-linear optimisation routine
%	through NETOPT. Note that for linear outputs, a single pass through
%	the  algorithm is all that is required, since the error function is
%	quadratic in the weights.  The error function value at the final set
%	of weights is returned in OPTIONS(8). Each row of X corresponds to
%	one input vector and each row of T corresponds to one target vector.
%
%	The optional parameters have the following interpretations.
%
%	OPTIONS(1) is set to 1 to display error values during training. If
%	OPTIONS(1) is set to 0, then only warning messages are displayed.  If
%	OPTIONS(1) is -1, then nothing is displayed.
%
%	OPTIONS(2) is a measure of the precision required for the value of
%	the weights W at the solution.
%
%	OPTIONS(3) is a measure of the precision required of the objective
%	function at the solution.  Both this and the previous condition must
%	be satisfied for termination.
%
%	OPTIONS(5) is set to 1 if an approximation to the Hessian (which
%	assumes that all outputs are independent) is used for softmax
%	outputs. With the default value of 0 the exact Hessian (which is more
%	expensive to compute) is used.
%
%	OPTIONS(14) is the maximum number of iterations for the IRLS
%	algorithm;  default 100.
%
%	See also
%	GLM, GLMERR, GLMGRAD
%

%	Copyright (c) Christopher M Bishop, Ian T Nabney (1996, 1997)

% Check arguments for consistency
errstring = consist(net, 'glm', x, t);
if ~errstring
  error(errstring);
end

if(~options(14))
  options(14) = 100;
end

display = options(1);

test = (options(2) | options(3));           % Do we need to test for termination?

ndata = size(x, 1);

inputs = [x ones(ndata, 1)];                % Add a column of ones for the bias 

% Use weighted iterative reweighted least squares (WIRLS)
e = ones(1, net.nin+1);
for n = 1:options(14)
    
    %switch net.actfn
      switch net.outfn
    case 'softmax'
      if n == 1        
        p = (t + (1/size(t, 2)))/2;    % Initialise model: ensure that row sum of p is one no matter
	    act = log(p./(1-p));           % how many classes there are
      end
      if options(5) == 1 | n == 1
        link_deriv = p.*(1-p);
        weights = sqrt(link_deriv);          % sqrt of weights
        if (min(min(weights)) < eps)
          fprintf(1, 'Warning: ill-conditioned weights in glmtrain\n')
          return
        end
        z = act + (t-p)./link_deriv;
        % Treat each output independently with relevant set of weights
        for j = 1:net.nout
          indep = inputs.*(weights(:,j)*e);
	      dep   = z(:,j).*weights(:,j);
	      temp  = indep\dep;
	      net.w1(:,j) = temp(1:net.nin);
	      net.b1(j)   = temp(net.nin+1);
        end
        [err, edata, eprior, p, act] = glmerr_weighted(net, x, t, eso_w);
        if n == 1
          errold = err;
          wold = netpak(net);
        else
          w = netpak(net);
        end
      else
	    % Exact method of calculation after w first initialised
	    % Start by working out Hessian
	    Hessian = glmhess_weighted(net, x, t, eso_w);
	    temp = p-t;
	    for m=1:ndata,
	       temp(m,:)=eso_w(m,1)*temp(m,:);
	    end
	    gw1 = x'*(temp);
	    gb1 = sum(temp, 1);
	    gradient = [gw1(:)', gb1];
	    % Now compute modification to weights
	    deltaw = -gradient*pinv(Hessian);
	    w = wold + alfa*deltaw;
	    net = glmunpak(net, w);
	    [err, edata, eprior, p] = glmerr_weighted(net, x, t, eso_w);
      end
    otherwise
      error(['Unknown activation function ', net.actfn]);
    end    % switch' end
    
   if options(1)==1
     fprintf(1, 'Cycle %4d Error %11.6f\n', n, err)
   end
   % Test for termination
   % Terminate if error increases
   if err >  errold
     errold = err;
     w = wold;
     options(8) = err;
     fprintf(1, 'Error has increased: terminating\n')
     return;
   end
   if test & n > 1
     if (max(abs(w - wold)) < options(2) & abs(err-errold) < options(3))
       options(8) = err;
       return;
     else
       errold = err;
       wold = w;
     end
   end
end

options(8) = err;
if (options(1) > 0)
  disp('Warning: Maximum number of iterations has been exceeded');
end
