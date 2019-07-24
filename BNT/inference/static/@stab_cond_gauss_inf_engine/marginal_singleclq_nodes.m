function marginal = marginal_singleclq_nodes(engine, i, query)
% MARGINAL_SINGLECLQ_NODES get the marginal distribution of nodes which is in a single clique
% marginal = marginal_singleclq_nodes(engine, i, query)

pot = struct(engine.clpot{i});
if isempty(pot.ctaildom)
    if i ~= engine.root
        p = parents(engine.jtree, i);
        tpot = direct_combine_pots(engine.clpot{i}, engine.seppot{p, i});
    else
        tpot = engine.clpot{i};
    end
    pot = marginalize_pot(tpot, query);
    
    marginal = pot_to_marginal(pot);
    marginal.T = normalise(marginal.T);
else
    [engine, clqtoroot] = push(engine, i, query);
    if clqtoroot == engine.root
        tpot = engine.clpot{clqtoroot};
    else
        p = parents(engine.jtree, clqtoroot);
        tpot = direct_combine_pots(engine.clpot{clqtoroot}, engine.seppot{p, clqtoroot});
    end
    pot = marginalize_pot(tpot, query);
    
    marginal = pot_to_marginal(pot);
    marginal.T = normalise(marginal.T);
end
                
