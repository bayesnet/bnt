function [ engine ] = jtree_inf_engine_rd( bnet, varargin )
% JTREE_INF_ENGINE2 Junction tree inference engine
% engine = jtree_inf_engine(bnet, ...)
%
% The following optional arguments can be specified in the form of name/value pairs:
% [default value in brackets]
%
% clusters  - a cell array of sets of nodes we want to ensure are in the same clique (in addition to families) [ {} ]
% root      - the root of the junction tree will be a clique that contains this set of nodes [N]
% stages    - stages{t} is a set of nodes we want to eliminate before stages{t+1}, ... [ {1:N} ]
%
% e.g., engine = jtree_inf_engine(bnet, 'maximize', 1);
%
% For more details on the junction tree algorithm, see
% - "Probabilistic networks and expert systems", Cowell, Dawid, Lauritzen and Spiegelhalter, Springer, 1999
% - "Inference in Belief Networks: A procedural guide", C. Huang and A. Darwiche, 
%      Intl. J. Approximate Reasoning, 15(3):225-263, 1996.
%
% modification of the calculus of cliques from JTREE_INF_ENGINE


% set default params
N = length( bnet.dag );

% Optional argument processing
[ b_verb, clusters, root, stages ] = jtree_inf_engine2_varargin_mgt( N, varargin );

% --- Verbose --- %
if b_verb
    fprintf( '/ ------------------------------- \\\n' );
    fprintf( '| jtree_inf_engine2 - verbose mode |\n' );
    time_i = datenum( clock );
end
% --------------- %

% Initialization
% ==============

% Class initialization
engine = init_fields;
engine = class( engine, 'jtree_inf_engine2', inf_engine( bnet ) );

% Default parameters
maximize = 0;
onodes = bnet.observed;

% Optional parameters given by user
% engine = set( engine, varargin{ : } );

% Building of the junction tree
% =============================

% Elimination ordering
% --------------------
% --- Verbose --- %
if b_verb
    fprintf( 'Compute elimination constraints ...' );
    time_i_cur = datenum( clock );
end
% --------------- %
porder = determine_elim_constraints( bnet, onodes );
strong = ~isempty( porder );
% --- Verbose --- %
if b_verb
    time_f_cur = datenum( clock );
    time_cur_str = datestr( time_f_cur - time_i_cur, 'HH:MM:SS' );
    fprintf( ' [Done] - elapsed time = %s\n', time_cur_str );
end
% --------------- %


% Moralization
% ------------
% --- Verbose --- %
if b_verb
    fprintf( 'Moralization ...' );
    time_i_cur = datenum( clock );
end
% --------------- %
ns = bnet.node_sizes( : );
ns( onodes ) = 1; % observed nodes have only 1 possible value
moral_graph = moralize( bnet.dag );
% --- Verbose --- %
if b_verb
    time_f_cur = datenum( clock );
    time_cur_str = datestr( time_f_cur - time_i_cur, 'HH:MM:SS' );
    fprintf( ' [Done] - elapsed time = %s\n', time_cur_str );
end
% --------------- %


% Junction tree building
% ----------------------
% --- Verbose --- %
if b_verb
    fprintf( 'Building the junction tree ...' );
    time_i_cur = datenum( clock );
end
% --------------- %
[ engine.jtree, root2, engine.cliques, B, w, elim_order ] = ...
    graph_to_jtree( moral_graph, ns, porder, stages, clusters );
% --- Verbose --- %
if b_verb
    time_f_cur = datenum( clock );
    time_cur_str = datestr( time_f_cur - time_i_cur, 'HH:MM:SS' );
    fprintf( ' [Done] - elapsed time = %s\n', time_cur_str );
end
% --------------- %


engine.cliques_bitv = B;
engine.clique_weight = w;
C = length( engine.cliques );
engine.clpot = cell(1,C);

% Separators computation
% ----------------------

% --- Verbose --- %
if b_verb
    fprintf( 'Separators computation ...' );
    time_i_cur = datenum( clock );
end
% --------------- %
% Compute the separators between connected cliques.
[ is, js ] = find( engine.jtree > 0 );
engine.separator = cell( C, C );
for k = 1:length( is )
  i = is( k ); j = js( k );
  % intersect(cliques{i}, cliques{j});
  engine.separator{ i, j } = find( B( i, : ) & B( j, : ) ); 
end
% --------------- %
if b_verb
    time_f_cur = datenum( clock );
    time_cur_str = datestr( time_f_cur - time_i_cur, 'HH:MM:SS' );
    fprintf( ' [Done] - elapsed time = %s\n', time_cur_str );
end
% --------------- %

% A node can be a member of many cliques, but is assigned to exactly one, to avoid
% double-counting its CPD. We assign node i to clique c if c is the "lightest" clique that
% contains i's family, so it can accomodate its CPD.

engine.clq_ass_to_node = zeros(1, N);
for i=1:N
  %c = clq_containing_nodes(engine, family(bnet.dag, i));
  % all selected columns must be 1
  clqs_containing_family = find( all( B( :, family( bnet.dag, i ) ), 2 ) );
  c = clqs_containing_family( ...
      argmin( w( clqs_containing_family ) ) );  
  engine.clq_ass_to_node( i ) = c; 
end

% Make the jtree rooted, so there is a fixed message passing order.
if strong
  % the last clique is guaranteed to be a strong root
  engine.root_clq = length( engine.cliques );
else
  % jtree_dbn_inf_engine requires the root to contain the interface.
  % This may conflict with the strong root requirement! *********** BUG *************
  engine.root_clq = clq_containing_nodes( engine, root );
  if engine.root_clq <= 0
    error( [ 'no clique contains ' num2str( root ) ] );
  end
end  

[ engine.jtree, engine.preorder, engine.postorder ] = ...
    mk_rooted_tree( engine.jtree, engine.root_clq );

% collect 
engine.postorder_parents = cell( 1, length(engine.postorder ) );
for n = engine.postorder( : )'
  engine.postorder_parents{ n } = parents( engine.jtree, n );
end
% distribute
engine.preorder_children = cell( 1, length( engine.preorder ) );
for n = engine.preorder( : )'
  engine.preorder_children{ n } = children( engine.jtree, n );
end

% --- Verbose --- %
if b_verb
    time_f = datenum( clock );
    time_str = datestr( time_f - time_i, 'HH:MM:SS' );
    fprintf( 'Elapsed time = %s\n', time_str' );
    fprintf( '\\ ------------------------------- /\n' );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ b_verb, clusters, root, stages ] = jtree_inf_engine2_varargin_mgt( N, parent_varargin )

% Number of variable arguments
nb_varargin = length( parent_varargin );

% Default parameters
b_verb = 0;
clusters = {};
root = N;
stages = { 1:N };

% Processing
for i = 1:2:nb_varargin
    
    arg_i = upper( parent_varargin{ i } );
    val_i = parent_varargin{ i + 1 };
    
    switch arg_i      
     case upper( 'EngineVerbose' )
      b_verb = val_i;
     case upper( 'Clusters' )
      clusters = val_i;
     case upper( 'Root' )
      root = val_i;
     case upper( 'Stages' )
      stages = val_i;
     otherwise
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function engine = init_fields()

engine.jtree = [];
engine.cliques = [];
engine.separator = [];
engine.cliques_bitv = [];
engine.clique_weight = [];
engine.clpot = [];
engine.clq_ass_to_node = [];
engine.root_clq = [];
engine.preorder = [];
engine.postorder = [];
engine.preorder_children = [];
engine.postorder_parents = [];
engine.maximize = [];
engine.evidence = [];

