function CPD = tabular_decision_node(sz, CPT)
% TABULAR_DECISION_NODE Represent the randomized policy over a discrete decision/action node as a table
% CPD = tabular_decision_node(sz, CPT)
%
% sz(1:end-1) is the sizes of the parents, sz(end) is the size of this node
% By default, CPT is set to the uniform random policy

if nargin==0
  % This occurs if we are trying to load an object from a file.
  CPD = init_fields;
  CPD = class(CPD, 'tabular_decision_node');
  return;
elseif isa(sz, 'tabular_decision_node')
  % This might occur if we are copying an object.
  CPD = sz;
  return;
end
CPD = init_fields;

if nargin < 2
  CPT = mk_stochastic(myones(sz)); 
else
  CPT = myreshape(CPT, sz);
end

CPD.CPT = CPT;
CPD.size = sz;

CPD = class(CPD, 'tabular_decision_node');

%%%%%%%%%%%

function CPD = init_fields()
% This ensures we define the fields in the same order 
% no matter whether we load an object from a file,
% or create it from scratch. (Matlab requires this.)

CPD.CPT = [];
CPD.size = [];
