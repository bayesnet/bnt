function vals = get_field(CPD, name)
% GET_PARAMS Get the parameters (fields) for a tabular_decision_node object
% vals = get_params(CPD, name)
%
% The following fields can be accessed
%
% policy - the table containing the policy
%
% e.g., policy = get_params(CPD, 'policy')

args = varargin;
nargs = length(args);
for i=1:2:nargs
  switch args{i},
   case 'policy',  vals =  CPD.CPT;
   otherwise,
    error(['invalid argument name ' args{i}]);
  end
end               
