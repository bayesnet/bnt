function engine = update_engine(engine, newCPDs)
% UPDATE_ENGINE Update the engine to take into account the new parameters (jtree_unrolled_dbn)
% engine = update_engine(engine, newCPDs)

engine.inf_engine = update_engine(engine.inf_engine, newCPDs);
engine.unrolled_engine = update_engine(engine.unrolled_engine, newCPDs);
                                                            
