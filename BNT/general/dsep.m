function sep = dsep(X, Y, S, G)
% DSEP Is X indep Y given S wrt DAG G?
% sep = dsep(X, Y, S, G)
%
% Instead of using the Bayes-Ball criterion, we see if S separates X and Y
% in the moralized ancestral graph.

conn = reachability_graph(G);
M = myunion(myunion(X, Y), S);
[A,junk] = find(conn(:, M));
A = unique(A);
A = myunion(A, M);
GM = moralize(G(A,A));
%sep = graph_separated(GM, X, Y, S);
sep = graph_separated(GM, find_equiv_posns(X,A), find_equiv_posns(Y,A), find_equiv_posns(S,A));
