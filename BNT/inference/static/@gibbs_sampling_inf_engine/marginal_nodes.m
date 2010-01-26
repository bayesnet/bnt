function [marginal, engine] = marginal_nodes(engine, nodes, varargin);
% MARGINAL_NODES Compute the marginal on the specified query nodes
% (gibbs_sampling_engine)
% [marginal, engine] = marginal_nodes(engine, nodes, ...)
%
% returns Pr(X(nodes) | X(observedNodes))
%
% The engine is also modified, and so it is returned as well, since
% Matlab doesn't support passing by reference(!)  So
% if you want to, for example, incrementally run gibbs for a few 100
% steps at a time, you should use the returned value.
%
% Optional arguments :
%
% 'reset_counts' is 1 if you want to reset the counts made in the
% past, and 0 otherwise (if the current query nodes are different
% from the previous query nodes, or if marginal_nodes has not been
% called before, reset_counts should be set to 1).
% By default it is 1.


reset_counts = 1;

if (nargin > 3)
  args = varargin;
  nargs = length(args);
  for i = 1:2:nargs
    switch args{i}
     case 'reset_counts'
      reset_counts = args{i+1};
     otherwise
      error(['Incorrect argument to gibbs_sampling_engine/' ...
	     ' marginal_nodes']);
    end
  end
end

% initialization stuff 
bnet = bnet_from_engine(engine);
slice_size = engine.slice_size;
hnodes = engine.hnodes;
onodes = engine.onodes;
nonqnodes = mysetdiff(1:slice_size, nodes);
gap = engine.gap;
burnin = engine.burnin;
T_max = engine.T;
ns = bnet.node_sizes(nodes);


% Cache the strides for the marginal table
marg_strides = [1 cumprod(ns(1:end-1))];
  
% Reset counts if necessary
if (reset_counts == 1) 
  %state = sample_bnet(bnet, 1, 0);
  %state = cell2num(sample_bnet(bnet, 'evidence', num2cell(engine.evidence)));
  state = cell2num(sample_bnet(bnet));
  state(onodes) = engine.evidence(onodes);
  if (length(ns) == 1)
    marginal_counts = zeros(ns(1),1);
  else
    marginal_counts = zeros(ns);
  end
  
% Otherwise, use the counts that have been stored in the engine  
else
  state = engine.state;
  state(onodes, :) = engine.evidence(onodes, :);
  marginal_counts = engine.marginal_counts;
end

if (engine.deterministic == 1)
  pos = 1;
  order = engine.order;
  orderSize = length(engine.order);
else
  sampling_dist = normalise(engine.sampling_dist);
end


for t = 1:(T_max*gap+burnin)

  % First, select node m to sample
  if (engine.deterministic == 1)
    m = engine.order(pos);
    pos = pos+1;
    if (pos > orderSize)
      pos = 1;
    end
  else
    m = my_sample_discrete(sampling_dist);
  end

  
  % If the node is observed, then don't bother resampling
  if (myismember(m, onodes))
    continue;
  end

  % Next, compute the posterior
  post = compute_posterior (bnet, state, m, engine.strides, engine.families, ...
			    engine.children, engine.CPT);
  state(m) = my_sample_discrete(post);

  % Now update our monte carlo estimate of the posterior
  % distribution on the query node 
  if ((mod(t-burnin, gap) == 0) & (t > burnin))

    vals = state(nodes);
    index = 1+marg_strides*(vals-1);
    marginal_counts(index) = marginal_counts(index)+1;
  end
end

% Store results for future computation.  Note that we store
% unnormalized counts
engine.state = state;
engine.marginal_counts = marginal_counts;

marginal.T = normalise(marginal_counts);


  
    











