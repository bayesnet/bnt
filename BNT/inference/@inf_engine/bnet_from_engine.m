function bnet = bnet_from_engine(engine)
% BNET_FROM_ENGINE Return the bnet structure stored inside the engine (inf_engine)
% bnet = bnet_from_engine(engine)

bnet = engine.bnet;

% We cannot write 'engine.bnet' without writing a 'subsref' function,
% since engine is an object with private parts.
% The bnet field should be the only thing external users of the engine should need access to.
% We do not pass bnet as a separate argument, since it could get out of synch with the one
% encoded inside the engine.
       
