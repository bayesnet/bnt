function val = get_params_kernel(K, name)
% GET_PARAMS_KERNEL Accessor function for a field (tabular_kernel)
% val = get_params_kernel(K, name)
%
% e.g., get_params_kernel(K, 'table')

switch name
 case 'table', val = K.table;
 otherwise,
  error(['invalid field name ' name]);
end
