function engine = enumerative_inf_engine(bnet)
% ENUMERATIVE_INF_ENGINE Inference engine for fully discrete BNs that uses exhaustive enumeration.
% engine = enumerative_inf_engine(bnet)


assert(isempty(bnet.cnodes));

% This is where we store stuff between enter_evidence and marginal_nodes
engine.evidence = [];

engine = class(engine, 'enumerative_inf_engine', inf_engine(bnet));
