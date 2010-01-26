function [pnodes, tnodes] = partition_dbn_nodes(intra, inter)
% PARTITION_DBN_NODES Divide the nodes into a DBN into persistent and transient.
% [pnodes, tnodes] = partition_dbn_nodes(intra, inter)
% Persistent nodes have children in the next time slice, transient nodes do not.
 
ss = length(intra);
pnodes = []; 
tnodes = [];
for i=1:ss
  cs = children(inter, i);
  if isempty(cs)
    tnodes = [tnodes i];
  else
    pnodes = [pnodes i];
  end
end
  
