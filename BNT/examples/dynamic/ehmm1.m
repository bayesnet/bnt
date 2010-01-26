% make the structure of an embedded HMM with 2 rows and 3 columns

% 1------------>2
% |\   \        | \  \
% 3->4->5       6->7->8

n = 8;
dag = zeros(n);
dag(1,[2 3 4 5])=1;
dag(2,[6 7 8])=1;
for i=3:4
  dag(i,i+1)=1;
end
for i=6:7
  dag(i,i+1)=1;
end
ns = 2*ones(1,n);
bnet = mk_bnet(dag,ns);
for i=1:n
  bnet.CPD{i}=tabular_CPD(bnet,i);
end
[jtree, root, cliques] =  graph_to_jtree(moralize(bnet.dag), ones(1,n), {}, {});
%[jtree, root, cliques, B, w, elim_order, moral_edges, fill_in_edges] = dag_to_jtree(bnet);

