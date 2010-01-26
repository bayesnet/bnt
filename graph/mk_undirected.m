function U = mk_undirected(G)

[nr nc] = size(G);
U = G;
for i=1:nr
  for j=1:nc
    if U(i,j)==1
      U(j,i) = 1;
    end
  end
end
