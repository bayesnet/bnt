function CPD = tabular_utility_node(bnet, node, T)
% TABULAR_UTILITY_NODE Represent a utility function as a table
% CPD = tabular_utility_node(bnet, node, T)
%
% node is the number of a node in this equivalence class.
% T is an optional argument (same shape as the CPT in tabular_CPD, but missing the last (child)
% dimension). By default, entries in T are chosen u.a.r. from 0:1 (using 'rand').

if nargin==0
  % This occurs if we are trying to load an object from a file.
  CPD = init_fields;
  clamp = 0;
  CPD = class(CPD, 'tabular_utility_node');
  return;
elseif isa(bnet, 'tabular_utility_node')
  % This might occur if we are copying an object.
  CPD = bnet;
  return;
end
CPD = init_fields;


ns = bnet.node_sizes;
ps = parents(bnet.dag, node);
sz = ns(ps);

if nargin < 3
  T = myrand(sz);
else
  T = myreshape(T, sz);
end

CPD.T = T;
CPD.sizes = sz;

CPD = class(CPD, 'tabular_utility_node');

%%%%%%%%%%%

function CPD = init_fields()
% This ensures we define the fields in the same order 
% no matter whether we load an object from a file,
% or create it from scratch. (Matlab requires this.)

CPD.T = [];
CPD.sizes = [];
