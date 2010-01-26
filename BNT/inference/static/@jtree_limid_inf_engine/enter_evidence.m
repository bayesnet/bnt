function engine = enter_evidence(engine, evidence, varargin)
% ENTER_EVIDENCE Add the specified evidence to the network (jtree_limid)
% engine = enter_evidence(engine, evidence, ...)
%
% evidence{i} = [] if if X(i) is hidden, and otherwise contains its observed value.
%
% The list below gives optional arguments [default value in brackets].      
%
% exclude - list of nodes whose potential will not be included in the joint [ [] ]
%
% e.g., engine = enter_evidence(engine, ev, 'exclude', 3)

exclude = [];

if nargin >= 3
  args = varargin;
  nargs = length(args);
  for i=1:2:nargs
    switch args{i},
     case 'exclude', exclude = args{i+1};
     otherwise,
      error(['invalid argument name ' args{i}]);
    end
  end
end
  
engine.exclude = exclude;
engine.evidence = evidence;
