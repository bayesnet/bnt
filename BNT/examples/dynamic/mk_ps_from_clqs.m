function mk_ps_from_clqs(dbn, T, cliques, dir)

% Draw multiple copies of the DBN,
% and indicate the nodes in each clique by shading the nodes.
% Generate a series of color postscript files,
% or, if dir=[], displays them to the screen and pauses.

if isempty(dir)
  print_to_file = 0;
else
  print_to_file = 1;
end

if print_to_file, cd(dir), end
flip = 1;
clf;
[dummyx, dummyy, h] = draw_dbn(dbn.intra, dbn.inter, flip, T, -1);

C = length(cliques);

% nodes = [];
% for i=1:C
%   cl = cliques{i};
%   nodes = [nodes cl(:)'];
% end
%nodes = unique(nodes);
ss = length(dbn.intra);
nodes = 1:(ss*T);

for c=1:C
  for i=cliques{c}
    set(h(i,2), 'facecolor', 'r'); 
  end
  rest = mysetdiff(nodes, cliques{c});
  for i=rest
    set(h(i,2), 'facecolor', 'w'); 
  end
  if print_to_file
    print(gcf, '-depsc', sprintf('clq%d.ps', c))
  else
   disp(['clique ' num2str(c) ' = ' num2str(cliques{c}) '; hit key for next'])
    pause
  end
end
