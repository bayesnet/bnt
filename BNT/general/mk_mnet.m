function mnet = mk_mnet(graph, node_sizes, cliques, potentials)
% MK_MNET Make a Markov network (Markov Random Field)
%
% mnet = mk_mnet(adj_mat, node_sizes, cliques, potentials)
%
% cliques{i} is a list of the nodes in clq i
% potentials{i} is a dpot object corresponding to the potential for clique i
%

mnet.markov_net = 1;
mnet.graph = graph;
mnet.node_sizes = node_sizes;
mnet.user_cliques = cliques;
mnet.user_potentials = potentials;
