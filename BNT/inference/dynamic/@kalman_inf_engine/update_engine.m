function engine = update_engine(engine, newCPDs)
% UPDATE_ENGINE Update the engine to take into account the new parameters (kalman)
% engine = update_engine(engine, newCPDs)

engine.inf_engine = update_engine(engine.inf_engine, newCPDs);
[engine.trans_mat, engine.trans_cov, engine.obs_mat, engine.obs_cov, engine.init_state, engine.init_cov] = ...
    dbn_to_lds(bnet_from_engine(engine));


