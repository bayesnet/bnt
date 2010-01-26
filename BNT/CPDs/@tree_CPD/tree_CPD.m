function CPD = tree_CPD(varargin)
%DTREE_CPD Make a conditional prob. distrib. which is a decision/regression tree.
%
% CPD =dtree_CPD() will create an empty tree.

if nargin==0
  % This occurs if we are trying to load an object from a file.
  CPD = init_fields;
  clamp = 0;
  CPD = class(CPD, 'tree_CPD', discrete_CPD(clamp, []));
  return;
elseif isa(varargin{1}, 'tree_CPD')
  % This might occur if we are copying an object.
  CPD = varargin{1};
  return;
end

CPD = init_fields;


clamped = 0;
fam_sz = [];
CPD = class(CPD, 'tree_CPD', discrete_CPD(clamped, fam_sz));


%%%%%%%%%%%

function CPD = init_fields()
% This ensures we define the fields in the same order 
% no matter whether we load an object from a file,
% or create it from scratch. (Matlab requires this.)

%init the decision tree set the root to null
CPD.tree.num_node = 0;
CPD.tree.root=1;
CPD.tree.nodes=[];

