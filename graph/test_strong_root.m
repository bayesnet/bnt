function strong = test_strong_root(jtree,cliques,dnodes,root)
% This function tests, whether root is a strong root of jtree. 
% The following parameters are used
% Input:
% jtree   An matrix with two colums. jtree(i,j) == jtree(j,i) is 1 if node 
%         i is connected with node j
% cliques Cells which contain the nodes in each clique
% dnodes  An array with the discrete nodes of the juntion tree.
% root    It is tested whether root is the strong root of the junction tree
% Output:
% strong  The output is 1 if root is the strong root of the junction tree jtree.
%         Please note, that the running intersection property is not tested. 
if isempty(dnodes)
     strong = 1;
     return;
end

children = find(jtree(root,:)==1);
i = 1;
strong = 1;
while (i <= length(children)) & (strong==1)
     child = children(i);
     jtree(child,root) = 0;
     jtree(root,child) = 0;
     sep = myintersect(cliques{child},cliques{root});
     diff = mysetdiff(cliques{child},cliques{root});
     if (mysubset(sep,dnodes) | isempty(myintersect(diff,dnodes)))
         strong = test_strong_root(jtree,cliques,dnodes,child);
     else
         strong = 0;
     end;
     i = i+1;
end

