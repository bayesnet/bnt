function bnet = fgraph_to_bnet(fg)
% FGRAPH_TO_BNET Convert a factor graph to a Bayes net
% bnet = fgraph_to_bnet(fg)
%
% We assume all factors are tabular_CPD.
% We create 1 dummy observed node for every factor.

N = fg.nvars + fg.nfactors;
vnodes = 1:fg.nvars;
fnodes = fg.nvars+1:N;
dag = zeros(N);
for x=1:fg.nvars
  dag(x, fnodes(fg.dep{x})) = 1;
end
ns = [fg.node_sizes ones(1, fg.nfactors)];
discrete = [fg.dnodes fnodes];
bnet = mk_bnet(dag, ns, 'discrete', discrete);
for x=1:fg.nvars
  bnet.CPD{x} = tabular_CPD(bnet, x, 'CPT', 'unif');
end
ev = cell(1, fg.nvars); % no evidence
for i=1:fg.nfactors
  f = fnodes(i);
  e = fg.equiv_class(i);
  pot = convert_to_pot(fg.factors{e}, 'd', fg.dom{i}, ev);
  m = pot_to_marginal(pot);
  bnet.CPD{f} = tabular_CPD(bnet, f, 'CPT', m.T);
end
  
  
