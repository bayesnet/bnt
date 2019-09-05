function smallpot = marginalize_pot(bigpot, keepdom, sumoverdom, nodesizes)
% MARGINALIZE_POT Marginalize a mpot onto a smaller domain.
% smallpot = marginalize_pot(bigpot, keep)

keepsize = sum(nodesizes(keepdom));
[A1, A2, B1, B2, C11, C12, C21, C22] = partition_matrix_vec_3(bigpot.A, bigpot.B, bigpot.C, keepdom, sumoverdom, nodesizes);
smallpot = scgcpot(keepsize, bigpot.ctailsize, bigpot.p, A1, B1, C11);