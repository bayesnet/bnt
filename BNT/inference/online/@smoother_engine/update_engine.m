function engine = update_engine(engine, newCPDs)
% UPDATE_ENGINE Update the engine to take into account the new parameters (smoother_engine).
% engine = update_engine(engine, newCPDs)

engine.tbn_engine = update_engine(engine.tbn_engine, newCPDs);
