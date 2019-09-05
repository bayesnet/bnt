%DEMGMM1 Demonstrate EM for Gaussian mixtures.
%
%	Description
%	This script demonstrates the use of the EM algorithm to fit a mixture
%	of Gaussians to a set of data using maximum likelihood. A colour
%	coding scheme is used to illustrate the evaluation of the posterior
%	probabilities in the E-step of the EM algorithm.
%
%	See also
%	DEMGMM2, DEMGMM3, DEMGMM4, GMM, GMMEM, GMMPOST
%

%	Copyright (c) Ian T Nabney (1996-2001)

mov = avifile('movies/gmm1.avi','fps',1 );

% Generate the data
randn('state', 0); rand('state', 0);
gmix = gmm(2, 2, 'spherical');
ndat1 = 20; ndat2 = 20; ndata = ndat1+ndat2;
gmix.centres =  [0.3 0.3; 0.7 0.7]; 
gmix.covars = [0.01 0.01];
x = gmmsamp(gmix, ndata);

h = figure;
hd = plot(x(:, 1), x(:, 2), '.g', 'markersize', 30);
hold on; axis([0 1 0 1]); axis square; set(gca, 'box', 'on');
ht = text(0.5, 1.05, 'Data', 'horizontalalignment', 'center');


% Set up mixture model
ncentres = 2; input_dim = 2;
mix = gmm(input_dim, ncentres, 'spherical');

% Initialise the mixture model
mix.centres = [0.2 0.8; 0.8, 0.2];
mix.covars = [0.01 0.01];

% Plot the initial model
ncirc = 30; theta = linspace(0, 2*pi, ncirc);
xs = cos(theta); ys = sin(theta);
xvals = mix.centres(:, 1)*ones(1,ncirc) + sqrt(mix.covars')*xs;
yvals = mix.centres(:, 2)*ones(1,ncirc) + sqrt(mix.covars')*ys;
hc(1)=line(xvals(1,:), yvals(1,:), 'color', 'r');
hc(2)=line(xvals(2,:), yvals(2,:), 'color', 'b');
set(ht, 'string', 'Initial Configuration');
figure(h);
mov = addframe(mov, getframe(gcf));
mov = addframe(mov, getframe(gcf));

% Initial E-step.
set(ht, 'string', 'E-step');
post = gmmpost(mix, x);
dcols = [post(:,1), zeros(ndata, 1), post(:,2)];
delete(hd); 
for i = 1 : ndata
  hd(i) = plot(x(i, 1), x(i, 2), 'color', dcols(i,:), ...
          'marker', '.', 'markersize', 30);
end

% M-step.
set(ht, 'string', 'M-step');
options = foptions; 
options(14) = 1; % A single iteration
options(1) = -1; % Switch off all messages, including warning
mix = gmmem(mix, x, options);
delete(hc);
xvals = mix.centres(:, 1)*ones(1,ncirc) + sqrt(mix.covars')*xs;
yvals = mix.centres(:, 2)*ones(1,ncirc) + sqrt(mix.covars')*ys;
hc(1)=line(xvals(1,:), yvals(1,:), 'color', 'r');
hc(2)=line(xvals(2,:), yvals(2,:), 'color', 'b');
figure(h);
mov = addframe(mov, getframe(gcf));
mov = addframe(mov, getframe(gcf));

% Loop over EM iterations.
numiters = 9;
for n = 1 : numiters

  set(ht, 'string', 'E-step');
  post = gmmpost(mix, x);
  dcols = [post(:,1), zeros(ndata, 1), post(:,2)];
  delete(hd); 
  for i = 1 : ndata
    hd(i) = plot(x(i, 1), x(i, 2), 'color', dcols(i,:), ...
                 'marker', '.', 'markersize', 30);
  end
  %pause(1)

  set(ht, 'string', 'M-step');
  [mix, options] = gmmem(mix, x, options);
  fprintf(1, 'Cycle %4d  Error %11.6f\n', n, options(8));
  delete(hc);
  xvals = mix.centres(:, 1)*ones(1,ncirc) + sqrt(mix.covars')*xs;
  yvals = mix.centres(:, 2)*ones(1,ncirc) + sqrt(mix.covars')*ys;
  hc(1)=line(xvals(1,:), yvals(1,:), 'color', 'r');
  hc(2)=line(xvals(2,:), yvals(2,:), 'color', 'b');
  pause(1)

  mov = addframe(mov, getframe(gcf));
end

mov = close(mov);
