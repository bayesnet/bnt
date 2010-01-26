function engine = push_pot_toclique(engine, clqtarget, clq, nodes)
% PUSH_POT push the variables in putshdom which is subset of clq to the target clique toword the root and get new engine
% engine = push_pot_toclique(engine, clqtarget, clq, nodes)
[engine, clqtoroot] = push_pot(engine, clq, nodes)
while clqtoroot ~= clqtarget
    [engine, clqtoroot] = push_pot(engine, clqtoroot, nodes)
end