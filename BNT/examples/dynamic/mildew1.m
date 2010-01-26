bnet = mk_mildew_dbn;

T = 4;
engine = {};
engine{end+1} = jtree_unrolled_dbn_inf_engine(bnet, T);
engine{end+1} = jtree_dbn_inf_engine(bnet); 
%engine{end+1} = hmm_inf_engine(bnet); % 8 is observed but has kids
engine{end+1} = smoother_engine(jtree_2TBN_inf_engine(bnet));
%engine{end+1} = smoother_engine(hmm_2TBN_inf_engine(bnet));

inf_time = cmp_inference_dbn(bnet, engine, T, 'check_ll', 0)
%learning_time = cmp_learning_dbn(bnet, engine, T)

S = struct(engine{1});
S1 = struct(S.unrolled_engine);
G = S1.jtree;
%graph_to_dot(G, 'directed', 0, 'leftright', 1, ...
%	     'filename', '/home/eecs/murphyk/WP/Thesis/Figures/Inf/Mildew/jtree.dot')
%!dot -Tps jtree.dot -o jtree.ps
% The resulting ps file cannot be converted using ps2pdf.

N = length(G);
for i=1:N
  for j=1:N
    if G(i,j)
      G(j,i)=1;
    end
  end
end
draw_graph(G)
