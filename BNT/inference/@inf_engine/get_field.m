function val = get_field(engine, name)
% GET_FIELD Get the value of a named field from a generic engine
% val = get_field(engine, name)
%
% The following fields can be accessed
%
% bnet
%
% e.g., bnet = get_field(engine, 'bnet')

switch name
 case 'bnet',      val = engine.bnet;
 otherwise,
  error(['invalid argument name ' name]);
end                                  
