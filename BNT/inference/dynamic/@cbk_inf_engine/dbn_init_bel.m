function engine = dbn_init_bel(engine)
% DBN_INIT_BEL Compute the initial belief state (bk)
% engine = dbn_init_bel(engine))

bnet = bnet_from_engine(engine);
ss = length(bnet.intra);
evidence = cell(1,ss);
engine = dbn_update_bel1(engine, evidence);
