function [engine, loglik] = enter_evidence(engine, evidence, varargin)
% ENTER_EVIDENCE Add the specified evidence to the network (jtree_online)
% [engine, loglik] = enter_evidence(engine, evidence, ...)
%
% evidence{i,t} = [] if if X(i,t) is hidden, and otherwise contains its observed value (scalar or column vector)
%
% The following optional arguments can be specified in the form of name/value pairs:
% [default value in brackets]
%
% maximize - if 1, does max-product instead of sum-product [0]
%

engine.maximize = 0;
args = varargin;
for i=1:2:length(args)
  switch args{i}
   case 'maximize', engine.maximize = args{i+1};
   otherwise, error(['unrecognized argument ' args{i}])
  end
end

[engine, loglik] = offline_smoother(engine, evidence);
