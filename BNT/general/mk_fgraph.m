function fg = mk_fgraph(G, node_sizes, factors, varargin)
% MK_FGRAPH Make a factor graph
% fg = mk_fgraph(G, node_sizes, factors, ...)
%
% A factor graph is a bipartite graph, with one side containing variables,
% and the other containing functions of (subsets of) these variables.
% For details, see "Factor Graphs and the Sum-Product Algorithm",
%  F. Kschischang and B. Frey and H-A. Loeliger,
%  IEEE Trans. Info. Theory, 2001
%
% G(i,j) = 1 if there is an arc from variable i to factor j
%
% node_sizes(i) is the number of values node i can take on,
%   or the length of node i if i is a continuous-valued vector.
%
% 'factors' is the list of factors (kernel functions)
%
% The list below gives optional arguments [default value in brackets].
% 
% equiv_class - equiv_class(i)=j  means factor node i gets its params from factors{j} [1:F]
% discrete - the list of nodes which are discrete random variables [1:N]
%
% e.g., fg = mk_fgraph(G, [2 2], {bnet.CPD{1},bnet.CPD{2}}, 'discrete', [1 2])

fg.G = G;
fg.node_sizes = node_sizes;
fg.factors = factors;
[fg.nvars fg.nfactors] = size(G);

% default values for parameters
fg.equiv_class = 1:fg.nfactors;
fg.dnodes = 1:fg.nvars;

if nargin >= 4
  args = varargin;
  nargs = length(args);
  for i=1:2:nargs
    switch args{i},
     case 'equiv_class', fg.equiv_class = args{i+1}; 
     case 'discrete',    fg.dnodes = args{i+1}; 
     otherwise,  
      error(['invalid argument name ' args{i}]);       
    end
  end
end

% so that determine_pot_type will work...
fg.utility_nodes = [];
%fg.decision_nodes = [];
%fg.chance_nodes = fg.nvars;

fg.dom = cell(1, fg.nfactors);
for f=1:fg.nfactors
  fg.dom{f} = find(G(:,f));
end
fg.dep = cell(1, fg.nvars);
for x=1:fg.nvars
  fg.dep{x} = find(G(x,:));
end
fg.cnodes = mysetdiff(1:fg.nvars, fg.dnodes);
