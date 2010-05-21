% =======
% display
% =======
%
% Description :
% -------------
%
% display function for class inf_engine.
%
% Syntax :
% --------
%
% [] = display( obj )
%
% Input(s) :
% ----------
%
% obj - class inf_engine
% An instance of the class inf_engine
%
% Output(s) :
% -----------
%
% Example(s) :
% ------------
%
% display( obj );
%
% Reference(s) :
% --------------
%
% See also :
% ----------
%
% disp
function [] = display( obj )

fprintf( '\n%s =\n\n', inputname( 1 ) ); 
disp( struct( obj ) );
    
% End of function

