function K = tabular_kernel(sz, table)
% TABULAR_KERNEL Make a table-based local kernel (discrete potential)
% K = tabular_kernel(sz, table)
%
% sz(i) is the number of values the i'th member of this kernel can have
% table is an optional array of size sz[1] x sz[2] x... [default: random]

if nargin==0
  % This occurs if we are trying to load an object from a file.
  K = init_fields;
  K = class(K, 'tabular_kernel');
  return;
elseif isa(sz, 'tabular_kernel')
  % This might occur if we are copying an object.
  K = sz;
  return;
end
K = init_fields;

if nargin < 2, table = myrand(sz); end

K.sz = sz;
K.table = table;

K = class(K, 'tabular_kernel');


%%%%%%%


function K = init_fields()
% This ensures we define the fields in the same order 
% no matter whether we load an object from a file,
% or create it from scratch. (Matlab requires this.)

K.sz = [];
K.table = [];



