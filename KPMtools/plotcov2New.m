% PLOTCOV2 - Plots a covariance ellipsoid with axes for a bivariate
%            Gaussian distribution.
%
% Usage:
%   [h, s] = plotcov2(mu, Sigma[, OPTIONS]);
% 
% Inputs:
%   mu    - a 2 x 1 vector giving the mean of the distribution.
%   Sigma - a 2 x 2 symmetric positive semi-definite matrix giving
%           the covariance of the distribution (or the zero matrix).
%
% Options:
%   'conf'      - a scalar between 0 and 1 giving the confidence
%                 interval (i.e., the fraction of probability mass to
%                 be enclosed by the ellipse); default is 0.9.
%   'num-pts'   - if the value supplied is n, then (n + 1)^2 points
%                 to be used to plot the ellipse; default is 20.
%   'label'     - if non-empty, a string that will label the
%                 ellipsoid (default: [])
%   'plot-axes' - a 0/1 flag indicating if the ellipsoid's axes
%                 should be plotted (default: 1)
%   'plot-opts' - a cell vector of arguments to be handed to PLOT3
%                 to contol the appearance of the axes, e.g., 
%                 {'Color', 'g', 'LineWidth', 1}; the default is {}
%   'fill-color' - a color specifier; is this is not [], the
%                  covariance ellipse is filled with this color
%                  (default: [])
% 
% Outputs:
%   h     - a vector of handles on the axis lines
%
% See also: PLOTCOV3

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

function [h, s] = plotcov2New(mu, Sigma, varargin)

h = [];
s = [];

if size(Sigma) ~= [2 2], error('Sigma must be a 2 by 2 matrix'); end
if length(mu) ~= 2, error('mu must be a 2 by 1 vector'); end

Sigma = checkpsd(Sigma);

[p, ...
 n, ...
 label, ...
 plot_axes, ...
 plot_opts, ...
 fill_color] = process_options(varargin, 'conf', 0.9, ...
			       'num-pts', 20, ...
			       'label', [], ...
			       'plot-axes', 1, ...
			       'plot-opts', {}, ...
			       'fill-color', []);
holding = ishold;
% Compute the Mahalanobis radius of the ellipsoid that encloses
% the desired probability mass.
k = conf2mahal(p, 2);
% Scale the covariance matrix so the confidence region has unit
% Mahalanobis distance.
Sigma = Sigma * k;
% The axes of the covariance ellipse are given by the eigenvectors of
% the covariance matrix.  Their lengths (for the ellipse with unit
% Mahalanobis radius) are given by the square roots of the
% corresponding eigenvalues.
[V, D] = eig(full(Sigma));
V = real(V);
D = real(D);
D = abs(D);

% Compute the points on the boundary of the ellipsoid.
t = linspace(0, 2*pi, n);
u = [cos(t(:))'; sin(t(:))'];
w = (V * sqrt(D)) * u;
z = repmat(mu(:), [1 n]) + w;
h = [h; plot(z(1, :), z(2, :), plot_opts{:})];
if (~isempty(fill_color))
  s = patch(z(1, :), z(2, :), fill_color);
end

% Plot the axes.
if (plot_axes)
  hold on;
  L = sqrt(diag(D));
  h = plot([mu(1); mu(1) + L(1) * V(1, 1)], ...
	   [mu(2); mu(2) + L(1) * V(2, 1)], plot_opts{:});
  h = [h; plot([mu(1); mu(1) + L(2) * V(1, 2)], ...
	       [mu(2); mu(2) + L(2) * V(2, 2)], plot_opts{:})];
end


if (~isempty(label))
  th = text(mu(1), mu(2), label);
  set(th, 'FontSize', 18);
  set(th, 'FontName', 'Times');
  set(th, 'FontWeight', 'bold');
  set(th, 'FontAngle', 'italic');
  set(th, 'HorizontalAlignment', 'center');
end

if (~holding & plot_axes) hold off; end
