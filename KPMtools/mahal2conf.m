% MAHAL2CONF - Translates a Mahalanobis distance into a confidence
%              interval.  Consider a multivariate Gaussian
%              distribution of the form
%
%   p(x) = 1/sqrt((2 * pi)^d * det(C)) * exp((-1/2) * MD(x, m, inv(C)))
%
%              where MD(x, m, P) is the Mahalanobis distance from x
%              to m under P:
%
%                 MD(x, m, P) = (x - m) * P * (x - m)'
%
%              A particular Mahalanobis distance k identifies an
%              ellipsoid centered at the mean of the distribution.
%              The confidence interval associated with this ellipsoid
%              is the probability mass enclosed by it.
%
%              If X is an d dimensional Gaussian-distributed vector,
%              then the Mahalanobis distance of X is distributed
%              according to the Chi-squared distribution with d
%              degrees of freedom.  Thus, the confidence interval is
%              determined by integrating the chi squared distribution
%              up to the Mahalanobis distance of the measurement.
%
% Usage:
% 
%   c = mahal2conf(m, d);
%
% Inputs:
%
%   m    - the Mahalanobis radius of the ellipsoid
%   d    - the number of dimensions of the Gaussian distribution
%
% Outputs:
%
%   c    - the confidence interval, i.e., the fraction of
%          probability mass enclosed by the ellipsoid with the
%          supplied Mahalanobis distance
%
% See also: CONF2MAHAL

% Copyright (C) 2002 Mark A. Paskin
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
% USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c = mahal2conf(m, d)

c = chi2cdf(m, d);