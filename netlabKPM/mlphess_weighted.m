function [h, hdata] = mlphess_weighted(net, x, t, eso_w, hdata)
%MLPHESS Evaluate the Hessian matrix for a multi-layer perceptron network.
%
%	Description
%	H = MLPHESS(NET, X, T) takes an MLP network data structure NET, a
%	matrix X of input values, and a matrix T of target values and returns
%	the full Hessian matrix H corresponding to the second derivatives of
%	the negative log posterior distribution, evaluated for the current
%	weight and bias values as defined by NET.
%
%	[H, HDATA] = MLPHESS(NET, X, T) returns both the Hessian matrix H and
%	the contribution HDATA arising from the data dependent term in the
%	Hessian.
%
%	H = MLPHESS(NET, X, T, HDATA) takes a network data structure NET, a
%	matrix X of input values, and a matrix T of  target values, together
%	with the contribution HDATA arising from the data dependent term in
%	the Hessian, and returns the full Hessian matrix H corresponding to
%	the second derivatives of the negative log posterior distribution.
%	This version saves computation time if HDATA has already been
%	evaluated for the current weight and bias values.
%
%	See also
%	MLP, HESSCHEK, MLPHDOTV, EVIDENCE
%

%	Copyright (c) Ian T Nabney (1996-9)

% Check arguments for consistency
errstring = consist(net, 'mlp', x, t);
if ~isempty(errstring);
  error(errstring);
end

if nargin == 4
  % Data term in Hessian needs to be computed
  hdata = datahess(net, x, t, eso_w);
end

[h, hdata] = hbayes(net, hdata);

% Sub-function to compute data part of Hessian
function hdata = datahess(net, x, t, eso_w)

hdata = zeros(net.nwts, net.nwts);

for v = eye(net.nwts);
  hdata(find(v),:) = mlphdotv_weighted(net, x, t, eso_w, v);
end

return
