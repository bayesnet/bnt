function CPD = tabular_decision_node(bnet, self, CPT)
% TABULAR_DECISION_NODE Represent a stochastic policy over a discrete decision/action node as a table
% CPD = tabular_decision_node(bnet, self, CPT)
%
% node is the number of a node in this equivalence class.
% CPT is an optional argument (see tabular_CPD for details); by default, it is the uniform policy.

if nargin==0
  % This occurs if we are trying to load an object from a file.
  CPD = init_fields;
  CPD = class(CPD, 'tabular_decision_node', discrete_CPD(1, []));
  return;
elseif isa(bnet, 'tabular_decision_node')
  % This might occur if we are copying an object.
  CPD = bnet;
  return;
end
CPD = init_fields;

ns = bnet.node_sizes;
fam = family(bnet.dag, self);
ps = parents(bnet.dag, self);
sz = ns(fam);

if nargin < 3
  CPT = mk_stochastic(myones(sz)); 
else
  CPT = myreshape(CPT, sz);
end

CPD.CPT = CPT;
CPD.sizes = sz; 

clamped = 1; % don't update using EM
CPD = class(CPD, 'tabular_decision_node', discrete_CPD(clamped, ns([ps self])));

%%%%%%%%%%%

function CPD = init_fields()
% This ensures we define the fields in the same order 
% no matter whether we load an object from a file,
% or create it from scratch. (Matlab requires this.)

CPD.CPT = [];
CPD.sizes = [];
