function fg = mk_fgraph_given_ev(G, node_sizes, factors, ev_CPD, evidence, varargin)
% MK_FGRAPH_GIVEN_EV Make a factor graph where each node has its own private evidence term
% fg = mk_fgraph(G, node_sizes, factors, ev_CPD, evidence, ...)
%
% G, node_sizes and factors are as in mk_fgraph, but they refer to the hidden nodes.
% ev_CPD{i} is a CPD for the i'th hidden node; this will be converted into a factor
% for node i using evidence{i}.
% We currently assume all hidden nodes are discrete, for simplicity.
%
% The list below gives optional arguments [default value in brackets].
% 
% equiv_class - equiv_class(i)=j  means factor node i gets its params from factors{j} [1:F]
% ev_equiv_class - ev_equiv_class(i)=j  means evidence node i gets its params from ev_CPD{j} [1:N]


N = length(node_sizes);
nfactors = length(factors);

% default values for parameters
eclass = 1:nfactors;
ev_eclass = 1:N;

if nargin >= 6
  args = varargin;
  nargs = length(args);
  for i=1:2:nargs
    switch args{i},
     case 'equiv_class', eclass = args{i+1}; 
     case 'ev_equiv_class', ev_eclass = args{i+1}; 
     otherwise,  
      error(['invalid argument name ' args{i}]);       
    end
  end
end

pot_type = 'd';
for x=1:N
  ev = cell(1,2); % cell 1 is the hidden parent, cell 2 is the observed child
  ev(2) = evidence(x); 
  dom = 1:2;
  F = convert_to_pot(ev_CPD{ev_eclass(x)}, pot_type, dom(:), ev);
  M = pot_to_marginal(F);
  %factors{end+1} = tabular_CPD('self', 1, 'ps', [], 'sz', node_sizes(x), 'CPT', M.T);
  factors{end+1} = mk_isolated_tabular_CPD(node_sizes(x), {'CPT', M.T});
end

E = max(eclass);
fg = mk_fgraph([G eye(N)], node_sizes, factors, 'equiv_class', [eclass E+1:E+N]);
