function [engine, loglik] = enter_evidence(engine, evidence, varargin)
% ENTER_EVIDENCE Add the specified evidence to the network (loopy_dbn)
% [engine, loglik] = enter_evidence(engine, evidence, ....)
%
% evidence{i,t} = [] if if X(i,t) is hidden, and otherwise contains its observed value (scalar or column vector)
%
% The following optional arguments can be specified in the form of name/value pairs:
% [default value in brackets]
%
% maximize - if 1, does max-product (not yet supported), else sum-product [0]
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
assert(~filter);

[engine.marginal, engine.msg, loglik] = enter_soft_ev(engine, evidence);
engine.evidence = evidence; % needed by marginal_nodes and marginal_family
