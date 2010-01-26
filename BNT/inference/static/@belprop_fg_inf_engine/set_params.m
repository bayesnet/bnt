function engine = set_params(engine, varargin)
% SET_PARAMS Set the parameters (fields) for a belprop_inf_engine object
% engine = set_params(engine, name/value pairs)
%
% The following optional arguments can be specified in the form of name/value pairs:
% e.g., engine = set_params(engine, 'tol', 1e-2, 'max_iter', 10)
%
% max_iter - max. num. loopy iterations 
% momentum - weight assigned to old message in convex combination 
% tol - tolerance used to assess convergence 
% maximize - 1 means use max-product, 0 means use sum-product

args = varargin{1};
nargs = length(args);
for i=1:2:nargs
  switch args{i},
   case 'max_iter', engine.max_iter = args{i+1};
   case 'momentum', engine.momentum = args{i+1};
   case 'tol',      engine.tol = args{i+1};
   case 'maximize', engine.maximize = args{i+1};
   otherwise,
    error(['invalid argument name ' args{i}]);
  end
end
