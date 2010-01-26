function [m, pot] = marginal_family(engine, query)
% MARGINAL_NODES Compute the marginal on the family of the specified node (jtree_limid)
% [m, pot] = marginal_family(engine, query)
%
% query should be a single decision node

bnet = bnet_from_engine(engine);
d = query;
assert(myismember(d, bnet.decision_nodes));
fam = family(bnet.dag, d);

clpot = init_clpot(bnet, engine.cliques, engine.clq_ass_to_node, engine.evidence, engine.exclude);

% collect to root (clique containing d) 
C = length(engine.cliques);
seppot = cell(C, C);    % separators are implicitely initialized to 1s
for n=engine.postorder{d}(1:end-1)
  for p=parents(engine.rooted_jtree{d}, n)
    %clpot{p} = divide_by_pot(clpot{n}, seppot{p,n}); % dividing by 1 is redundant
    seppot{p,n} = marginalize_pot(clpot{n}, engine.separator{p,n});
    clpot{p} = multiply_by_pot(clpot{p}, seppot{p,n});
  end
end

root = engine.clq_ass_to_node(d);
assert(root == engine.postorder{d}(end));
pot = marginalize_pot(clpot{root}, fam);
m = pot_to_marginal(pot);

%%%%%%%%%%%


function clpot = init_clpot(bnet, cliques, clq_ass_to_node, evidence, exclude)

% Set the clique potentials to all 1s
C = length(cliques);
clpot = cell(1, C);
ns = bnet.node_sizes;
for i=1:C
  clpot{i} = upot(cliques{i}, ns(cliques{i}));
end

N = length(bnet.dag);
nodes = mysetdiff(1:N, exclude);

for n=nodes(:)'
  fam = family(bnet.dag, n);
  e = bnet.equiv_class(n);
  c = clq_ass_to_node(n);
  pot = convert_to_pot(bnet.CPD{e}, 'u', fam(:), evidence);
  clpot{c} = multiply_by_pot(clpot{c}, pot);
end
