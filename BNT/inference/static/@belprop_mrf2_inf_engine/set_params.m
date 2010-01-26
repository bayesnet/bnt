function engine = set_params(engine, varargin)
% SET_PARAMS Modify parameters of the inference engine
% engine = set_params(engine, 'param1',val1, 'param2',val2, ...)
%
% Parameter names are listed below.
%
% max_iter - max. num. iterations 
% momentum - weight assigned to old message in convex combination
%            (useful for damping oscillations) 
% tol      - tolerance used to assess convergence
% verbose - 1 means print error at every iteration [0]

[engine.max_iter, engine.momentum, engine.tol, engine.verbose] = ...
    process_options('max_iter', engine.max_iter, 'momentum', engine.momentum, ...
		    'tol', engine.tol, 'verbose', engine.verbose);
