function mrf2 = mk_mrf2(adj_mat, pot)
% MK_MRF2 Make a Markov random field with pairwise potentials
% function mrf2 = mk_mrf2(adj_mat, pot)
%
% pot{i,j}(k1,k2)

mrf2.adj_mat = adj_mat;
mrf2.pot = pot;
