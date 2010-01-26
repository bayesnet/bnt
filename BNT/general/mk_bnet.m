function bnet = mk_bnet(dag, node_sizes, varargin)
% MK_BNET Make a Bayesian network.
%
% BNET = MK_BNET(DAG, NODE_SIZES, ...) makes a graphical model with an arc from i to j iff DAG(i,j) = 1.
% Thus DAG is the adjacency matrix for a directed acyclic graph.
% The nodes are assumed to be in topological order. Use TOPOLOGICAL_SORT if necessary.
%
% node_sizes(i) is the number of values node i can take on,
%   or the length of node i if i is a continuous-valued vector.
% node_sizes(i) = 1 if i is a utility node.
% 
% Below are the names of optional arguments [and their default value in brackets].
% Pass as 'PropertyName1', PropertyValue1, 'PropertyName2', PropertyValue2, ...
% 
% discrete - the list of nodes which are discrete random variables [1:N]
% equiv_class - equiv_class(i)=j  means node i gets its params from CPD{j} [1:N]
% observed - the list of nodes which will definitely be observed in every case [ [] ]
% 'names' - a cell array of strings to be associated with nodes 1:n [{}]
%    This creates an associative array, so you write e.g.
%     'evidence(bnet.names{'bar'}) = 42' instead of  'evidence(2} = 42' 
%     assuming names = { 'foo', 'bar', ...}.
%
% e.g., bnet = mk_bnet(dag, ns, 'discrete', [1 3])
%
% For backwards compatibility with BNT2, you can also specify the parameters in the following order
%   bnet = mk_bnet(dag, node_sizes, discrete_nodes, equiv_class)

n = length(dag);

% default values for parameters
bnet.equiv_class = 1:n;
bnet.dnodes = 1:n; % discrete 
bnet.observed = [];
bnet.names = {};

if nargin >= 3
  args = varargin;
  nargs = length(args);
  if ~isstr(args{1})
    if nargs >= 1, bnet.dnodes = args{1}; end
    if nargs >= 2, bnet.equiv_class = args{2}; end
  else    
    for i=1:2:nargs
      switch args{i},
       case 'equiv_class', bnet.equiv_class = args{i+1}; 
       case 'discrete',    bnet.dnodes = args{i+1}; 
       case 'observed',    bnet.observed = args{i+1}; 
       case 'names',  bnet.names = assocarray(args{i+1}, num2cell(1:n)); 
       otherwise,  
	error(['invalid argument name ' args{i}]);       
      end
    end
  end
end

bnet.observed = sort(bnet.observed); % for comparing sets
bnet.hidden = mysetdiff(1:n, bnet.observed(:)');
bnet.hidden_bitv = zeros(1,n);
bnet.hidden_bitv(bnet.hidden) = 1;
bnet.dag = dag;
bnet.node_sizes = node_sizes(:)';

bnet.cnodes = mysetdiff(1:n, bnet.dnodes);
% too many functions refer to cnodes to rename it to cts_nodes - 
% We hope it won't be confused with chance nodes!

bnet.parents = cell(1,n);
for i=1:n
  bnet.parents{i} = parents(dag, i);
end

E = max(bnet.equiv_class);
mem = cell(1,E);
for i=1:n
  e = bnet.equiv_class(i);
  mem{e} = [mem{e} i];
end
bnet.members_of_equiv_class = mem;

bnet.CPD = cell(1, E);

bnet.rep_of_eclass = zeros(1,E);
for e=1:E
  mems = bnet.members_of_equiv_class{e};
  bnet.rep_of_eclass(e) = mems(1);
end

directed = 1;
if ~acyclic(dag,directed)
  error('graph must be acyclic')
end

bnet.order = topological_sort(bnet.dag);
