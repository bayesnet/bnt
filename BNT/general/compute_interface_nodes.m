function [interface, persist, transient] = compute_interface_nodes(intra, inter)
% COMPUTE_INTERFACE_NODES Find the nodes in a DBN that represent a sufficient statistic
% [interface, persist, transient] = compute_interface_nodes(intra, inter)
%
% The interface nodes are all those that has an incoming temporal arc,
% or which are parents of such nodes.
% If the parents are in the previous slice, this just means they have an
% outgoing temporal arc.
% (The parents of nodes with incoming temporal arcs are needed
% because moralization will bring them into the clique.)
%
% The persisent nodes are all those that have one or more incoming temporal arc.
% The transient nodes are all the non-persistent.
%
% See U. Kjaerulff, "dHugin: A computational system for dynamic
% time-sliced Bayesian networks", Intl. J. Forecasting (11) 89-111, 1995

n = length(intra);
interface = [];
persist = [];
% any nodes with incoming arcs
for u=1:n
  if any(inter(:,u))
    interface = [interface u];
    persist = [persist u];
  end
end
% Any nodes which are parents of nodes with incoming arcs
for u=1:n
  cs = children(intra, u);
  if any(inter(:, cs))
    interface = [interface u];
  end
  %cs = children(inter, u);
  % if ~isempty(myintersect(cs, persist))
  %  interface = [interface u];
  %end
end
interface = unique(interface);
persist = unique(persist);
transient = mysetdiff(1:n, persist);

