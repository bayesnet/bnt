function mpe = find_mpe(engine, evidence, varargin)
% FIND_MPE Find the most probable explanation of the data  (belprop_fg)
% function mpe = find_mpe(engine, evidence,...)
%
% evidence{i} = [] if X(i) is hidden, and otherwise contains its observed value (scalar or column vector).
%
% This finds the marginally most likely value for each hidden node,
% and may give the wrong results even if the graph is acyclic,
% unless you set break_ties = 1.
%
% The following optional arguments can be specified in the form of name/value pairs:
% [default value in brackets]
%
% break_ties is optional. If 1, we will force ties to be broken consistently
%  by calling enter_evidence N times. (see Jensen96, p106) Default = 1.

break_ties = 1;

% parse optional params
args = varargin;
nargs = length(args);
for i=1:2:nargs
  switch args{i},
   case 'break_ties',    break_ties = args{i+1}; 
   otherwise,  
    error(['invalid argument name ' args{i}]);       
  end
end

engine = enter_evidence(engine, evidence, 'maximize', 1);

observed = ~isemptycell(evidence);
evidence = evidence(:); % hack to handle unrolled DBNs
N = length(evidence);
mpe = cell(1,N);
for i=1:N
  m = marginal_nodes(engine, i);
  % observed nodes are all set to 1 inside the inference engine, so we must undo this
  if observed(i)
    mpe{i} = evidence{i};
  else
    mpe{i} = argmax(m.T);
    if break_ties
      evidence{i} = mpe{i};                             
      [engine, ll] = enter_evidence(engine, evidence, 'maximize', 1);  
    end
  end
end

