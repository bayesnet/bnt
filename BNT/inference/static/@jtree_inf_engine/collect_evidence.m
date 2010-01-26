function [clpot, seppot] = collect_evidence(engine, clpot, seppot)
% COLLECT_EVIDENCE Do message passing from leaves to root (children then parents)
% [clpot, seppot] = collect_evidence(engine, clpot, seppot)

for n=engine.postorder %postorder(1:end-1)
  for p=engine.postorder_parents{n}
    %clpot{p} = divide_by_pot(clpot{n}, seppot{p,n}); % dividing by 1 is redundant
    seppot{p,n} = marginalize_pot(clpot{n}, engine.separator{p,n}, engine.maximize);
    clpot{p} = multiply_by_pot(clpot{p}, seppot{p,n});
  end
end

