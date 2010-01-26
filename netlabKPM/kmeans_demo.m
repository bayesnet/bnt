function kmeans_demo()

% Generate T points from K=5 1D clusters, and try to recover the cluster
% centers using k-means.
% Requires BNT, netlab and the matlab stats toolbox v4.

K = 5;
ndim = 1;
true_centers = 1:K;
sigma = 1e-6;
T = 100;
% data(t,:) is the t'th data point
data = zeros(T, ndim); 
% ndx(t) = i means the t'th data point is sample from cluster i
%ndx = sample_discrete(normalise(ones(1,K)));
ndx = [1*ones(1,20) 2*ones(1,20) 3*ones(1,20) 4*ones(1,20) 5*ones(1,20)];
for t=1:T
  data(t) = sample_gaussian(true_centers(ndx(t)), sigma, 1);
end
plot(1:T, data, 'x')



% set the centers randomly from Gauss(0)
mix = gmm(ndim, K, 'spherical');
h = plot_centers_as_lines(mix, [], T);

if 0
% Place initial centers at K data points chosen at random, but add some noise
choose_ndx = randperm(T);
choose_ndx = choose_ndx(1:K);
init_centers = data(choose_ndx) + sample_gaussian(0, 0.1, K);
mix.centres = init_centers;
h = plot_centers_as_lines(mix, h, T);
end

if 0
% update centers using netlab k-means
options = foptions;
niter = 10;
options(14) = niter;
mix = gmminit(mix, data, options);
h = plot_centers_as_lines(mix, h, T);
end

% use matlab stats toolbox k-means with multiple restarts
nrestarts = 5;
[idx, centers] = kmeans(data, K, 'replicates', nrestarts, ...
			'emptyAction', 'singleton', 'display', 'iter');
mix.centres = centers;
h = plot_centers_as_lines(mix, h, T);

% fine tune with EM; compute covariances of each cluster
options = foptions;
niter = 20;
options(1) = 1; % display cost fn at each iter
options(14) = niter;
mix = gmmem(mix, data, options);
h = plot_centers_as_lines(mix, h, T);

%%%%%%%%%
function h = plot_centers_as_lines(mix, h, T)

K = mix.ncentres;
hold on
if isempty(h)
  for k=1:K
    h(k)=line([0 T], [mix.centres(k) mix.centres(k)]);
  end
else
  for k=1:K
    set(h(k), 'xdata', [0 T], 'ydata', [mix.centres(k) mix.centres(k)]);
  end
end
hold off

