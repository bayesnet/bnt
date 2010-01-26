function bnet = bnet_from_engine(engine)
% BNET_FROM_ENGINE Return the bnet structure stored inside the engine (smoother_engine)
% bnet = bnet_from_engine(engine)

bnet = bnet_from_engine(engine.tbn_engine);
