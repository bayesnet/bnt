function marginal = dbn_marginal_from_bel(engine, i)
% DBN_MARGINAL_FROM_BEL Compute the marginal on a node given the current belief state (bk_ff_hmm)
% marginal = dbn_marginal_from_bel(engine, i)

marginal = pot_to_marginal(engine.bel_marginals{i});
