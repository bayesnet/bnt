function K = tabular_kernel(fg, self)
% TABULAR_KERNEL Make a table-based local kernel (discrete potential)
% K = tabular_kernel(fg, self)
%
% fg is a factor graph
% self is the number of a representative domain
%
% Use 'set_params_kernel' to adjust the following fields
%   table - a q[1]xq[2]x... array, where q[i] is the number of values for i'th node
%       in this domain [default: random values from [0,1], which need not sum to 1]


if nargin==0
  % This occurs if we are trying to load an object from a file.
  K = init_fields;
  K = class(K, 'tabular_kernel');
  return;
elseif isa(fg, 'tabular_kernel')
  % This might occur if we are copying an object.
  K = fg;
  return;
end
K = init_fields;

ns = fg.node_sizes;
dom = fg.doms{self};
% we don't store the actual domain since it may vary due to parameter tieing
K.sz = ns(dom);
K.table = myrand(K.sz);

K = class(K, 'tabular_kernel');


%%%%%%%


function K = init_fields()
% This ensures we define the fields in the same order 
% no matter whether we load an object from a file,
% or create it from scratch. (Matlab requires this.)

K.table = [];
K.sz = [];


