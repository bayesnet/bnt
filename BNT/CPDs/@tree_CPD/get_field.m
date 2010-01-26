function val = get_params(CPD, name)
% GET_PARAMS Get the parameters (fields) for a tabular_CPD object
% val = get_params(CPD, name)
%
% The following fields can be accessed
%
% cpt       - the CPT
%
% e.g., CPT = get_params(CPD, 'cpt')

switch name
 case 'cpt',      val = CPD.CPT;
 case 'tree',     val = CPD.tree;
 otherwise,
  error(['invalid argument name ' name]);
end
