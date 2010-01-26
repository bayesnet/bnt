function CPD = root_CPD(bnet, self, val)
% ROOT_CPD Make a conditional prob. distrib. which has no parameters.
% CPD = ROOT_CPD(BNET, NODE_NUM, VAL)
%
% The node must not have any parents and is assumed to always be observed.
% It is a way of modelling exogenous inputs to a model.
% VAL is the value to which the root is clamped (default: [])


if nargin==0
  % This occurs if we are trying to load an object from a file.
  CPD = init_fields;
  CPD = class(CPD, 'root_CPD', generic_CPD(1));
  return;
elseif isa(bnet, 'root_CPD')
  % This might occur if we are copying an object.
  CPD = bnet;
  return;
end
CPD = init_fields;


if nargin < 3, val = []; end

ns = bnet.node_sizes;
ps = parents(bnet.dag, self);
if ~isempty(ps)
  error('root CPDs should have no parents')
end

CPD.self = self;
CPD.val = val;
CPD.sizes = ns(self);

clamped = 1;
CPD = class(CPD, 'root_CPD', generic_CPD(clamped));


%%%%%%%%%%%

function CPD = init_fields()
% This ensures we define the fields in the same order 
% no matter whether we load an object from a file,
% or create it from scratch. (Matlab requires this.)

CPD.self = [];
CPD.val = [];
CPD.sizes = [];
