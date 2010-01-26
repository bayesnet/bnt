function [A, names] = mk_adj_mat(connections, names, topological)
% MK_ADJ_MAT Make a directed adjacency matrix from a list of connections between named nodes.
%
% A = mk_adj_mat(connections, name)
% This is best explaine by an example:
%   names = {'WetGrass', 'Sprinkler', 'Cloudy', 'Rain'}; 
%   connections = {'Cloudy', 'Sprinkler'; 'Cloudy', 'Rain'; 'Sprinkler', 'WetGrass'; 'Rain', 'WetGrass'}; 
% adds the arcs C -> S, C -> R, S -> W, R -> W. Node 1 is W, 2 is S, 3 is C, 4 is R.
%
% [A, names] = mk_adj_mat(connections, name, 1)
% The last argument of 1 indicates that we should topologically sort the nodes (parents before children).
% In the example, the numbering becomes: node 1 is C, 2 is R, 3 is S, 4 is W
% and the return value of names gets permuted to {'Cloudy', 'Rain', 'Sprinkler', 'WetGrass'}.
% Note that topological sorting the graph is only possible if it has no directed cycles.

if nargin < 3, topological = 0; end
  
n=length(names);
A=zeros(n);
[nr nc] = size(connections);
for r=1:nr
  from = strmatch(connections{r,1}, names, 'exact');
  assert(~isempty(from));
  to = strmatch(connections{r,2}, names, 'exact');
  assert(~isempty(to));
  %fprintf(1, 'from %s %d to %s %d\n', connections{r,1}, from, connections{r,2}, to);
  A(from,to) = 1;
end

if topological
  order = topological_sort(A); 
  A = A(order, order); 
  names = names(order); 
end


