function sep = graph_separated(G, X, Y, S)

G2 = G;
G2(S,:) = 0;
G2(:,S) = 0;
conn = reachability_graph(G2);
conn2 = conn(X,Y);
sep = all(conn2(:)==0);

