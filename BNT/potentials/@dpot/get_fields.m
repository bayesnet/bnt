function val = get_params(pot, name)
% GET_PARAMS Accessor function for a field (dpot)
% val = get_params(pot, name)
%
% e.g., get_params(pot, 'table') or 'domain'

switch name
 case 'table', val = pot.T;
 case 'domain', val = pot.domain;
 otherwise,
  error(['invalid field name ' name]);
end

