% ====
% disp
% ====
%
% Description :
% -------------
%
% disp function for class inf_engine.
%
% Syntax :
% --------
%
% [] = disp( obj )
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
% disp( obj );
%
% Reference(s) :
% --------------
%
% See also :
% ----------
%
% display
function [] = disp( obj )

disp( struct( obj ) );
    
% End of function
