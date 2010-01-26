function gdl = mk_gdl_graph(G, domains, node_sizes, kernels, varargin)
% MK_GDL_GRAPH Make a GDL (generalized distributed law) graph
% gdl = mk_gdl_graph(G, domains, node_sizes, kernels, ...)
%
% A GDL graph is like a moralized, but untriangulated, Bayes net:
% each "node" represents a domain with a corresponding kernel function.
% For details, see "The Generalized Distributive Law", Aji and McEliece,
% IEEE Trans. Info. Theory, 46(2): 325--343, 2000
% 
% G(i,j) = 1 if there is an (undirected) edge between domains i,j
%
% domains{i} is the domain of node i
%
% node_sizes(i) is the number of values node i can take on,
%   or the length of node i if i is a continuous-valued vector.
% node_sizes(i) = 1 if i is a utility node.
%
% kernels is the list of kernel functions
%
% The list below gives optional arguments [default value in brackets].
% 
% equiv_class - equiv_class(i)=j  means factor node i gets its params from factors{j} [1:F]
% discrete - the list of nodes which are discrete random variables [1:N]
% chance   - the list of nodes which are random variables [1:N]
% decision - the list of nodes which are decision nodes [ [] ]
% utility  - the list of nodes which are utility nodes [ [] ]


ns = node_sizes;
N = length(domains);
vars = [];
for i=1:N
  vars = myunion(vars, domains{i});
end
Nvars  = length(vars);

gdl.equiv_class = 1:length(kernels);
gdl.chance_nodes = 1:Nvars;
gdl.utility_nodes = [];
gdl.decision_nodes = [];
gdl.dnodes = 1:Nvars;

if nargin >= 5
  args = varargin;
  nargs = length(args);
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


gdl.G = G;
gdl.vars = vars;
gdl.doms = domains;
gdl.node_sizes = node_sizes;
gdl.cnodes = mysetdiff(vars, gdl.dnodes);
gdl.kernels = kernels;
gdl.type = 'gdl';

% Compute a bit vector representation of the set of domains
% dom_bitv(i,j) = 1 iff variable j occurs in domain i
gdl.dom_bitv = zeros(N, length(vars));
for i=1:N
  gdl.dom_bitv(i, domains{i}) = 1;
end

% compute the interesection of the domains on either side of each edge (separating set)
gdl.sepset = cell(N, N);
gdl.nbrs = cell(1,N);
for i=1:N
  nbrs = neighbors(G, i);
  gdl.nbrs{i} = nbrs;
  for j = nbrs(:)'
    gdl.sepset{i,j} = myintersect(domains{i}, domains{j});
  end
end


