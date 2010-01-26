function [int, persist, transient] = compute_interface_nodes(intra, inter)
% COMPUTE_INTERFACE_NODES Find the nodes in a DBN that represent a sufficient statistic
% [int, persist, transient] = compute_interface_nodes(intra, inter)
%
% The interface nodes are all those that has an incoming temporal arc,
% or which have a child which has an incoming temporal arc,
% where a temporal arc means one coming from the previous slice.
% (The parents of nodes with incoming temporal arcs are needed
% because moralization will bring them into the clique.)
%
% The persisent nodes are all those that have one or more incoming temporal arc.
% The transient nodes are all the non-persistent.
%
% See U. Kjaerulff, "dHugin: A computational system for dynamic
% time-sliced Bayesian networks", Intl. J. Forecasting (11) 89-111, 1995

n = length(intra);
int = [];
persist = [];
for u=1:n
  if any(inter(:,u))
    int = [int u];
    persist = [persist u];
  end
  if any(inter(:, children(intra, u)))
    int = [int u];
  end
end
int = unique(int);
persist = unique(persist);
transient = mysetdiff(1:n, persist);
