function bnet = mutilate_bnet(bnet, nodes, vals)
% MUTILATE_BNET Clamp nodes to specific values (perform a surgical intervention)
% bnet = mutilate_bnet(bnet, nodes, vals)
%
% We make all the clamped nodes roots.

ns = bnet.node_sizes;
for i=1:length(nodes)
  X = nodes(i);
  x = vals(i);
  bnet.dag(:,X) = 0;
  bnet.CPD{X} = root_CPD(bnet, X, x);
end
