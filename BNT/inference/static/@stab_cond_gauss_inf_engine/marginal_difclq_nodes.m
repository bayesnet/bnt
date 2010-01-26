function marginal = marginal_difclq_nodes(engine, query_nodes)
% MARGINAL_DIFCLQ_NODES get the marginal distribution of nodes which is not in a single clique
% marginal = marginal_difclq_nodes(engine, query_nodes)

keyboard
num_clique = length(engine.cliques);
B = engine.cliques_bitv;
clqs_containnodes = [];
for i=1:length(query_nodes)
    node = query_nodes(i);
    tnodes = find(all(B(:, node), 2));
    clqs_containnodes = myunion(clqs_containnodes, tnodes);
end
% get all cliques contains query nodes

% get the minimal sub tree in junction which contains these cliques and the node closest to the root of jtree
[subtree, nroot_node] = min_subtree_conti_nodes(engine.jtree, engine.root, clqs_containnodes);
if ~mysubset(query_nodes, engine.cliques{nroot_node});
    % if query nodes is not all memers of the clique closest to the root clique performe push operation
    engine = push_tree(engine, subtree, query_nodes, nroot_node);
end

if ~(nroot_node == engine.root)
    % if the clique closest to the root clique is not the root clique we must direct combine the 
    % potential with the potential stored in separator toward to root
    p = parents(engine.jtree, nroot_node);
    tpot = direct_combine_pots(engine.clpot{nroot_node}, engine.seppot{p, nroot_node});
else
    tpot = engine.clpot{nroot_node};
end

pot = marginalize_pot(tpot, query_nodes);
marginal = pot_to_marginal(pot);
marginal.T = normalise(marginal.T);



function engine = push_tree(engine, tree, query_nodes, inode)
% PUSH_TREE recursive perform push opeartion on tree
% engine = push_tree(engine, tree, query_nodes, inode)

cs = children(tree, inode);
for i = 1:length(cs)
    node = cs(i);
    push_tree(engine, tree, query_nodes, node);
    push_dom = myintersect(engine.cliques{node}, query_nodes);
    [engine, clqtoroot] = push(engine, node, push_dom);
end







