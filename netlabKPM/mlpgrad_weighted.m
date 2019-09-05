function [g, gdata, gprior] = mlpgrad_weighted(net, x, t, eso_w)
%MLPGRAD Evaluate gradient of error function for 2-layer network.
%
%	Description
%	G = MLPGRAD(NET, X, T) takes a network data structure NET  together
%	with a matrix X of input vectors and a matrix T of target vectors,
%	and evaluates the gradient G of the error function with respect to
%	the network weights. The error funcion corresponds to the choice of
%	output unit activation function. Each row of X corresponds to one
%	input vector and each row of T corresponds to one target vector.
%
%	[G, GDATA, GPRIOR] = MLPGRAD(NET, X, T) also returns separately  the
%	data and prior contributions to the gradient. In the case of multiple
%	groups in the prior, GPRIOR is a matrix with a row for each group and
%	a column for each weight parameter.
%
%	See also
%	MLP, MLPPAK, MLPUNPAK, MLPFWD, MLPERR, MLPBKP
%

%	Copyright (c) Ian T Nabney (1996-9)

% Check arguments for consistency
errstring = consist(net, 'mlp', x, t);
if ~isempty(errstring);
  error(errstring);
end
[y, z] = mlpfwd(net, x);
temp = y - t;
ndata = size(x, 1);
for m=1:ndata,
      delout(m,:)=eso_w(m,1)*temp(m,:);
end
clear temp;
gdata = mlpbkp(net, x, z, delout);

[g, gdata, gprior] = gbayes(net, gdata);
