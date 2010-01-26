function engine = belprop_inf_engine(bnet, varargin) 
% BELPROP_INF_ENGINE Make a loopy belief propagation inference engine
% engine = belprop_inf_engine(bnet, ...)
%
% This is like pearl_inf_engine, except it uses potential objects,
% instead of lambda/pi structs. Hence it is slower.
%
% The following optional arguments can be specified in the form of name/value pairs:
% [default in brackets]
%
% protocol - 'tree' means send messages up then down the tree,
%            'parallel' means use synchronous updates ['parallel']
% max_iter - max. num. iterations [ 2*num_nodes ]
% momentum - weight assigned to old message in convex combination (useful for damping oscillations) [0]
% tol      - tolerance used to assess convergence [1e-3]
% maximize - 1 means use max-product, 0 means use sum-product [0]
% filename - name of file to write beliefs to after each iteration within enter_evidence [ [] ]
%
% e.g., engine = belprop_inf_engine(bnet, 'maximize', 1, 'max_iter', 10)

% gdl = general distributive law
engine.gdl = bnet_to_gdl(bnet);

% set default params
N = length(engine.gdl.G);
engine.protocol = 'parallel';
engine.max_iter = 2*N;
engine.momentum = 0;
engine.tol = 1e-3;
engine.maximize = 0;
engine.filename = [];
engine.fid = [];

args = varargin;
nargs = length(args);
for i=1:2:nargs
  switch args{i},
   case 'max_iter', engine.max_iter = args{i+1};
   case 'momentum', engine.momentum = args{i+1};
   case 'tol',      engine.tol = args{i+1};
   case 'protocol', engine.protocol = args{i+1};
   case 'filename', engine.filename = args{i+1};
   otherwise,
    error(['invalid argument name ' args{i}]);
  end
end


if strcmp(engine.protocol, 'tree')
  % Make a rooted tree, so there is a fixed message passing order.
  root = N;
  [engine.tree, engine.preorder, engine.postorder, height, cyclic] = mk_rooted_tree(engine.gdl.G, root);
  assert(~cyclic);
end

% store results computed by enter_evidence here
engine.marginal_domains = cell(1, N);

engine.niter = [];

engine = class(engine, 'belprop_inf_engine', inf_engine(bnet));

%%%%%%%%%

function gdl = bnet_to_gdl(bnet)

gdl.G = mk_undirected(bnet.dag);
N = length(bnet.dag);
gdl.doms = cell(1,N);
for i=1:N
  gdl.doms{i} = family(bnet.dag, i);
end 

% Compute a bit vector representation of the set of domains
% dom_bitv(i,j) = 1 iff variable j occurs in domain i
gdl.dom_bitv = zeros(N, N);
for i=1:N
  gdl.dom_bitv(i, gdl.doms{i}) = 1;
end
   
% compute the interesection of the domains on either side of each edge (separating set)
gdl.sepset = cell(N, N);
gdl.nbrs = cell(1,N);
for i=1:N
  nbrs = neighbors(gdl.G, i);
  gdl.nbrs{i} = nbrs;
  for j = nbrs(:)'
    gdl.sepset{i,j} = myintersect(gdl.doms{i}, gdl.doms{j});
  end
end  
