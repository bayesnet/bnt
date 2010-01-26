function bnet15 = mk_slice_and_half_dbn(bnet, int)
% function bnet = mk_slice_and_half_dbn(bnet, int)
% function bnet = mk_slice_and_half_dbn(bnet, int)
%
% Create a "1.5 slice" jtree, containing the interface nodes of slice 1
% and all the nodes of slice 2
% To keep the node numbering the same, we simply disconnect the non-interface nodes
% from slice 1, and set their size to 1.
% We do this to speed things up, and so that the likelihood is computed correctly.
% We do not need to do
% this if we just want to compute marginals (i.e., we can include nodes whose potentials will
% be left as all 1s).

intra15 = bnet.intra;
ss = length(bnet.intra);
nonint = mysetdiff(1:ss, int);
for i=nonint(:)'
  intra15(:,i) = 0;
  intra15(i,:) = 0;
  %assert(~any(bnet.inter(i,:)))
end
dag15 = [intra15      bnet.inter;
	 zeros(ss)    bnet.intra];
ns = bnet.node_sizes(:);
ns(nonint) = 1; % disconnected nodes get size 1
obs_nodes = [bnet.observed(:) bnet.observed(:)+ss];
bnet15 = mk_bnet(dag15, ns, 'discrete', bnet.dnodes, 'equiv_class', bnet.equiv_class(:), ...
		 'observed', obs_nodes(:));
