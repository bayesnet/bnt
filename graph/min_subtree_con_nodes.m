function [subtree, nroot_node] = min_subtree_con_nodes(jtree, root, nodes)
%min_subtree_con_nodes get the minimum subtree of tree which contains the nodes

if isempty(jtree) | isempty(nodes)
    subtree = [];
    nroot_node = [];
    return;
end

rnodes = min_subtree_nodes(jtree, nodes);
nea_node = nearest_node(jtree, root, nodes);
node_num = length(jtree);
subtree = zeros(node_num);
subtree(rnodes, rnodes) = jtree(rnodes, rnodes);
nroot_node = nea_node;


function rnodes = min_subtree_nodes(tree, nodes)
rnodes = [];
if isempty(tree) | isempty(nodes)
    return
end

rnodes = nodes(1);
newnodes = neighbors(tree, nodes(1));
while ~mysubset(nodes, rnodes)
    swapnodes = newnodes;
    newnodes = [];
    added = 0;
    for i=1:length(swapnodes)
        inode = swapnodes(i);
        tnodes = myunion(inode, rnodes);
        if mysubset(nodes, tnodes)
            added = 1;
            break;
        end
        nns = neighbors(tree, inode);
        add_nodes = mysetdiff(nns, tnodes);
        newnodes = myunion(newnodes, add_nodes);
    end
    if added
        rnodes = tnodes;
    else
        rnodes = myunion(rnodes, newnodes);
    end
end

function nea_node = nearest_node(tree, inode, nodes)
if myismember(inode, nodes)
    nea_node = inode;
    return;
end
cs = children(tree, inode);
for i = 1:length(cs)
    n = cs(i);
    nea_node = nearest_node(tree, n, nodes);
end
    


