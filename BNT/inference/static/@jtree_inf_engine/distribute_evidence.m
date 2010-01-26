function [clpot, seppot] = distribute_evidence(engine, clpot, seppot)
% DISTRIBUTE_EVIDENCE Do message passing from root to leaves (parents then children)
% [clpot, seppot] = distribute_evidence(engine, clpot, seppot)

for n=engine.preorder
  for c=engine.preorder_children{n}
    clpot{c} = divide_by_pot(clpot{c}, seppot{n,c}); 
    seppot{n,c} = marginalize_pot(clpot{n}, engine.separator{n,c}, engine.maximize);
    clpot{c} = multiply_by_pot(clpot{c}, seppot{n,c});
  end
end
