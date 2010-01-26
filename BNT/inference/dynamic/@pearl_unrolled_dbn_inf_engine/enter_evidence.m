function [engine, loglik, niter] = enter_evidence(engine, evidence, varargin)
% ENTER_EVIDENCE Add the specified evidence to the network (loopy_unrolled_dbn)
% [engine, loglik, niter] = enter_evidence(engine, evidence, ....)
%
% evidence{i,t} = [] if if X(i,t) is hidden, and otherwise contains its observed value (scalar or column vector)
%
% The following optional arguments can be specified in the form of name/value pairs:
% [default value in brackets]
%
% maximize - if 1, does max-product (not yet supported), else sum-product [0]
% filename - as in loopy_pearl
%
% e.g., engine = enter_evidence(engine, ev, 'maximize', 1)

maximize = 0;
filename = engine.filename;

if nargin >= 2
  args = varargin;
  nargs = length(args);
  for i=1:2:nargs
    switch args{i},
     case 'maximize', maximize = args{i+1};
     case 'filename', filename = args{i+1};
    end
  end
end


[ss T] = size(evidence);
if T ~= engine.T
  bnetT = dbn_to_bnet(bnet_from_engine(engine), T);
  engine.unrolled_engine = pearl_inf_engine(bnetT, 'protocol', engine.protocol, ...
					    'max_iter', engine.max_iter_per_slice * T, ...
					    'tol', engine.tol, 'momentum', engine.momentum);
  engine.T = T;
end
[engine.unrolled_engine, loglik, niter] = enter_evidence(engine.unrolled_engine, evidence(:), ...
						  'maximize', maximize, 'filename', filename);


