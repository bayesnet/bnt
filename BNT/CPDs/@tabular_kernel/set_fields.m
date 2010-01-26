function K = set_params_kernel(K, name, val)
% SET_PARAMS_KERNEL Accessor function for a field (table_kernel)
% K = set_params_kernel(K, name, val)
%
% e.g., K = set_params_kernel(K, 'table', rand(2,3,2)) for a kernel on 3 nodes with 2,3,2 values each

% We should check if the arguments are valid...

switch name
 case 'table', K.table = val;
 otherwise,
  error(['invalid field name ' name]);
end
