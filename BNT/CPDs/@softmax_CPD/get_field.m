function val = get_params(CPD, name)
% GET_PARAMS Get the parameters (fields) for a softmax_CPD object
% val = get_params(CPD, name)
%
% The following fields can be accessed
%
% weights - W(X,Y,Q)
% offset  - b(Y,Q)
%
% e.g., W = get_params(CPD, 'weights')

[W, b] = extract_params(CPD);
switch name
 case 'weights',   val = W;
 case 'offset',    val = b;
 otherwise,
  error(['invalid argument name ' name]);
end                
