% ===
% set
% ===
%
% Description :
% -------------
%
% Writing method for the class' attributes. 
% 
% Syntax :
% --------
%
% obj = set( obj, attrib_name_1, val_1, ..., attrib_name_n, val_n );
%
% Input(s) :
% ----------
%
% obj - class inf_engine
% An instance of the class inf_engine
%
% > Optionals : 'attrib_name'/value form
% For details about the attribute's names, see the constructor
% (inf_engine.m)
%
% Output(s) :
% -----------
%
% obj - class inf_engine
% The updated object
%
% Example(s) :
% ------------
%
% obj = set( obj, 'attrib_1', val_1, ..., 'attrib_n', val_n );
%
% Reference(s) :
% --------------
%
% See also :
% ----------
%
% inf_engine
% get
function [ engine ] = set( engine, varargin )

engine.inf_engine = set( engine.inf_engine, varargin{ : } );

% Processing
% Number of variable arguments
nb_varargin = length( varargin );

for i = 1:2:nb_varargin
    
    arg_i = upper( varargin{ i } );
    val_i = varargin{ i + 1 };
    
    switch arg_i      
     case upper( 'maximize' )
      engine.maximize = val_i;
     otherwise
%       warning_str = [ 'inf_engine - set : attribute ' ...
% 		      arg_i ' doesn''t exist' ];
		      
%       warning( warning_str );
    end
    
end
    
% End of function