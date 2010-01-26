function CPD = set_params(CPD, varargin)
% SET_PARAMS Set the parameters (fields) for a tabular_decision_node object
% CPD = set_params(CPD, name/value pairs)
%
% The following optional arguments can be specified in the form of name/value pairs:
%
% policy - the table containing the policy
%
% e.g., CPD = set_params(CPD, 'policy', T)

args = varargin;
nargs = length(args);
for i=1:2:nargs
  switch args{i},
   case 'policy',   CPD.CPT = args{i+1};
   otherwise,
    error(['invalid argument name ' args{i}]);
  end
end               
