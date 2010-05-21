function b = isdag(G)
% b = isdag(G)
% 
b = sum(sum(G.*G'));        % How many undirected arcs ? (x2)
b=~b & ~isempty(G);
if b
  M = expm(full(G)) - eye(length(G)); M = (M>0);
  b = isempty(find(sum(sum(eye(length(G)).*M)))); % is there no cycle ?
end
