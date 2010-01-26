function bnet = mk_limid(dag, node_sizes, varargin)
% MK_LIMID Make a limited information influence diagram
%
% BNET = MK_LIMID(DAG, NODE_SIZES, ...) 
% DAG is the adjacency matrix for a directed acyclic graph.
% The nodes are assumed to be in topological order. Use TOPOLOGICAL_SORT if necessary.
% For decision nodes, the parents must explicitely include all nodes
% on which it can depends, in contrast to the implicit no-forgetting assumption of influence diagrams.
% (For details, see "Representing and solving decision problems with limited information",
%    Lauritzen and Nilsson, Management Science, 2001.)
%
% node_sizes(i) is the number of values node i can take on,
%   or the length of node i if i is a continuous-valued vector.
% node_sizes(i) = 1 if i is a utility node.
% 
% The list below gives optional arguments [default value in brackets].
% 
% chance   - the list of nodes which are random variables [1:N]
% decision - the list of nodes which are decision nodes [ [] ]
% utility  - the list of nodes which are utility nodes [ [] ]
% equiv_class - equiv_class(i)=j  means node i gets its params from CPD{j} [1:N]
%
% e.g., limid = mk_limid(dag, ns, 'chance', [1 3], 'utility', [2])

n = length(dag);

% default values for parameters
bnet.chance_nodes = 1:n;
bnet.equiv_class = 1:n;
bnet.utility_nodes = [];
bnet.decision_nodes = [];
bnet.dnodes = 1:n; % discrete 

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
       case 'chance',      bnet.chance_nodes = args{i+1}; 
       case 'utility',     bnet.utility_nodes = args{i+1}; 
       case 'decision',    bnet.decision_nodes = args{i+1}; 
       case 'discrete',    bnet.dnodes = args{i+1}; 
        otherwise,  
	error(['invalid argument name ' args{i}]);       
      end
    end
  end
end
 
bnet.limid = 1;

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

% for e=1:E
%   i = bnet.members_of_equiv_class{e}(1); % pick arbitrary member
%   switch type{e}
%     case 'tabular',  bnet.CPD{e} = tabular_CPD(bnet, i);
%     case 'gaussian', bnet.CPD{e} = gaussian_CPD(bnet, i);
%     otherwise, error(['unrecognized CPD type ' type{e}]);
%   end
% end

directed = 1;
if ~acyclic(dag,directed)
  error('graph must be acyclic')
end

bnet.order = topological_sort(bnet.dag);
