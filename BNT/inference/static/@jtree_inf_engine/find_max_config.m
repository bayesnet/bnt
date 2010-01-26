function [mpe, clpot, seppot] = find_max_config(engine, clpot, seppot, evidence)
% FIND_MAX_CONFIG Backwards pass of Viterbi fro jtree
% function [mpe, clpot, seppot] = find_max_config(engine, clpot, seppot, evidence)
% See Cowell99 p98

bnet = bnet_from_engine(engine);
nnodes = length(bnet.dag);
mpe = cell(1, nnodes);
maximize = 1;

c = engine.root_clq;
pot = struct(clpot{c}); % violate object privacy
dom = pot.domain;
[indices, clpot{c}] = find_most_prob_entry(clpot{c});
mpe(dom) = num2cell(indices);

for n=engine.preorder
  for c=engine.preorder_children{n}
    clpot{c} = divide_by_pot(clpot{c}, seppot{n,c}); 
    seppot{n,c} = marginalize_pot(clpot{n}, engine.separator{n,c}, maximize);
    clpot{c} = multiply_by_pot(clpot{c}, seppot{n,c});
    
    pot = struct(clpot{c}); % violate object privacy
    dom = pot.domain;
    [indices, clpot{c}] = find_most_prob_entry(clpot{c});
    mpe(dom) = num2cell(indices);
  end
end

obs_nodes = find(~isemptycell(evidence));
% indices for observed nodes will be 1 - need to overwrite these
mpe(obs_nodes) = evidence(obs_nodes);



