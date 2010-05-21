% ===
% get
% ===
%
% Description :
% -------------
%
% Reading method for the class' attributes. 
% 
% Syntax :
% --------
%
% val = get( obj, attribute_name );
%
% Input(s) :
% ----------
%
% obj - class inf_engine
% An instance of the class inf_engine
%
% attribute_name - string
% The name of a class' attribute
% For details about the attribute's names, see the constructor
% (inf_engine.m)
%
% Output(s) :
% -----------
%
% val - any
% The value of the specified attribute
%
% Example(s) :
% ------------
%
% val = get( obj, 'attribute_name' );
%
% Reference(s) :
% --------------
%
% See also :
% ----------
%
% inf_engine
% set
function [ val ] = get( engine, attribute_name )

val = [];

% Processing 
val = get( engine.inf_engine, attribute_name );

% Processing    
arg = upper( attribute_name );

switch arg      
 case upper( 'evidence' )
  val = engine.evidence;
 case upper( 'jtree' )
  val = engine.jtree;
 case upper( 'cliques' )
  val = engine.cliques;
  case upper( 'separator' )
  val = engine.separator;
  case upper( 'cliques_bitv' )
  val = engine.cliques_bitv;
  case upper( 'clique_weight' )
  val = engine.clique_weight;
  case upper( 'clpot' )
  val = engine.clpot;
  case upper( 'clq_ass_to_node' )
  val = engine.clq_ass_to_node;
 case upper( 'root_clq' )
  val = engine.root_clq;
  case upper( 'preorder' )
  val = engine.preorder;
  case upper( 'postorder' )
  val = engine.postorder;
  case upper( 'preorder_children' )
  val = engine.preorder_children;
  case upper( 'postorder_parents' )
  val = engine.postorder_parents;
  case upper( 'maximize' )
  val = engine.maximize;
  case upper( 'evidence' )
  val = engine.evidence;
             
 otherwise
  %   warning_str = [ 'inf_engine - get : attribute ' ...
  % 		  arg ' doesn''t exist' ];
  
  %   warning( warning_str );
end

% End of function

