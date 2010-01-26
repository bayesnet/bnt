function [trans_mat, trans_cov, obs_mat, obs_cov, init_state, init_cov] = dbn_to_lds(bnet)
% DBN_TO_LDS Compute the Linear Dynamical System parameters from the Gaussian DBN.
% [trans_mat, trans_cov, obs_mat, obs_cov, init_state, init_cov] = dbn_to_lds(bnet)

onodes = bnet.observed;
ss = length(bnet.intra);
num_nodes = ss*2;
assert(isequal(bnet.cnodes_slice, 1:ss));
[W,D,mu] = extract_params_from_gbn(bnet);

hnodes = mysetdiff(1:ss, onodes);
bs = bnet.node_sizes(:); % block sizes

obs_mat = W(block(hnodes,bs), block(onodes,bs))';
u = block(onodes,bs);
obs_cov = D(u,u);

trans_mat = W(block(hnodes,bs), block(hnodes + ss, bs))';
u = block(hnodes + ss, bs);
trans_cov = D(u,u);

u = block(hnodes,bs);
init_cov = D(u,u);
init_state = mu(u);


