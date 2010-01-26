function engine = update_engine(engine, newCPDs)
% UPDATE_ENGINE Update the engine to take into account the new parameters (hmm)
% engine = update_engine(engine, newCPDs)

%engine.inf_engine.bnet.CPD = newCPDs;
engine.inf_engine = update_engine(engine.inf_engine, newCPDs);
[engine.startprob, engine.transprob, engine.obsprob] = dbn_to_hmm(bnet_from_engine(engine));

