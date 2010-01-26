
%bnet = mk_uffe_dbn;
bnet = mk_mildew_dbn;
ss = length(bnet.intra);

% construct jtree from 1.5 slice DBN

int = compute_fwd_interface(bnet.intra, bnet.inter);
bnet15 = mk_slice_and_half_dbn(bnet, int);

% use unconstrained elimination,
% but force there to be a clique containing both interfaces
clusters = {int, int+ss};
jtree_engine = jtree_inf_engine(bnet15, 'clusters', clusters, 'root', int+ss);
S=struct(jtree_engine)
in_clq = clq_containing_nodes(jtree_engine, int);
out_clq = clq_containing_nodes(jtree_engine, int+ss)


% Also make a jtree from slice 1
bnet1 = mk_bnet(bnet.intra1, bnet.node_sizes_slice);
jtree_engine1 = jtree_inf_engine(bnet1, 'clusters', {int}, 'root', int);
S1=struct(jtree_engine1)
