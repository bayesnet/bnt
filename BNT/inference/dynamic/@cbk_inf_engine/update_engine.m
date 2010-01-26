function engine = update_engine(engine, newCPDs)
% UPDATE_ENGINE Update the engine to take into account the new parameters (bk)
% engine = update_engine(engine, newCPDs)

engine.inf_engine = update_engine(engine.inf_engine, newCPDs);
engine.sub_engine = update_engine(engine.sub_engine, newCPDs);

bnet = bnet_from_engine(engine);
eclass1 = bnet.equiv_class(:,1);
engine.sub_engine1 = update_engine(engine.sub_engine1, newCPDs(1:max(eclass1)));

