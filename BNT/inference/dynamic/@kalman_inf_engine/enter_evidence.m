function [engine, loglik] = enter_evidence(engine, evidence, varargin)
% ENTER_EVIDENCE Add the specified evidence to the network (kalman)
% [engine, loglik] = enter_evidence(engine, evidence, ...)
%
% evidence{i,t} = [] if if X(i,t) is hidden, and otherwise contains its observed value (scalar or column vector)
%
% The following optional arguments can be specified in the form of name/value pairs:
% [default value in brackets]
%
% maximize - if 1, does max-product (same as sum-product for Gaussians!), else sum-product [0]
% filter -   if 1, do filtering, else smoothing [0]
%
% e.g., engine = enter_evidence(engine, ev, 'maximize', 1)

maximize = 0;
filter = 0;

% parse optional params
args = varargin;
nargs = length(args);
if nargs > 0
  for i=1:2:nargs
    switch args{i},
     case 'maximize', maximize = args{i+1}; 
     case 'filter', filter = args{i+1}; 
     otherwise,  
      error(['invalid argument name ' args{i}]);       
    end
  end
end

assert(~maximize);

bnet = bnet_from_engine(engine);
n = length(bnet.intra);
onodes = bnet.observed;
hnodes = mysetdiff(1:n, onodes);
T = size(evidence, 2);
ns = bnet.node_sizes;
O = sum(ns(onodes));
data = reshape(cat(1, evidence{onodes,:}), [O T]);

A = engine.trans_mat;
C = engine.obs_mat;
Q = engine.trans_cov;
R = engine.obs_cov;
init_x = engine.init_state;
init_V = engine.init_cov;

if filter
  [x, V, VV, loglik] = kalman_filter(data, A, C, Q, R, init_x, init_V);
else
  [x, V, VV, loglik] = kalman_smoother(data, A, C, Q, R, init_x, init_V);
end

  
% Wrap the posterior inside a potential, so it can be marginalized easily
engine.one_slice_marginal = cell(1,T);
engine.two_slice_marginal = cell(1,T);
ns(onodes) = 0;
ns(onodes+n) = 0;
ss = length(bnet.intra);
for t=1:T
  dom = (1:n);
  engine.one_slice_marginal{t} = mpot(dom+(t-1)*ss, ns(dom), 1, x(:,t), V(:,:,t));
end
% for t=1:T-1
%   dom = (1:(2*n));
%   mu = [x(:,t); x(:,t)];
%   Sigma = [V(:,:,t) VV(:,:,t+1)';
% 	   VV(:,:,t+1) V(:,:,t+1)];
%   engine.two_slice_marginal{t} = mpot(dom+(t-1)*ss, ns(dom), 1, mu, Sigma);
% end
for t=2:T
  %dom = (1:(2*n));
  current_slice = hnodes;
  next_slice = hnodes + ss;
  dom = [current_slice next_slice];   
  mu = [x(:,t-1); x(:,t)];
  Sigma = [V(:,:,t-1) VV(:,:,t)';
	   VV(:,:,t) V(:,:,t)];
  engine.two_slice_marginal{t-1} = mpot(dom+(t-2)*ss, ns(dom), 1, mu, Sigma);
end
