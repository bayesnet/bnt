function engine = dbn_init_bel(engine)
% DBN_INIT_BEL Compute the initial belief state (bk_ff_hmm)
% engine = dbn_init_bel(engine)

engine.bel = engine.prior(:);
