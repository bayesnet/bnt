function engine = belprop_mrf2_inf_engine(mrf2, varargin) 
% BELPROP_MRF2_INF_ENGINE Belief propagation for MRFs with discrete pairwise potentials
% engine = belprop_mrf2_inf_engine(mrf2, ...)
%
% This is like belprop_inf_engine, except it is designed for mrf2, so is much faster.
%
% [ ... ] = belprop_mrf2_inf_engine(..., 'param1',val1, 'param2',val2, ...)
% allows you to specify optional parameters as name/value pairs.
% Parameters modifying behavior of enter_evidence are below [default value in brackets]
%
% max_iter - max. num. iterations [ 5*nnodes]
% momentum - weight assigned to old message in convex combination
%            (useful for damping oscillations) [0]
% tol      - tolerance used to assess convergence [1e-3]
% verbose - 1 means print error at every iteration [0]
%
% Parameters can be changed later using set_params 


% The advantages of pairwise potentials are
% (1) we can compute messages using vector-matrix multiplication
% (2) we can easily specify the parameters: one potential per edge
% In contrast, potentials on larger cliques are more complicated to deal with.


nnodes = length(mrf2.adj_mat);

[engine.max_iter, engine.momentum, engine.tol, engine.verbose] = ...
    process_options(varargin, 'max_iter', [], 'momentum', 0, 'tol', 1e-3, ...
		   'verbose', 0);

if isempty(engine.max_iter) % no user supplied value, so compute default
  engine.max_iter = 5*nnodes;
  %if acyclic(mrf2.adj_mat, 0) --- can be very slow!
  %  engine.max_iter = nnodes;
  %else
  %  engine.max_iter = 5*nnodes;
  %end
end

engine.bel = cell(1, nnodes); % store results of enter_evidence here
engine.mrf2 = mrf2;

engine = class(engine, 'belprop_mrf2_inf_engine');


