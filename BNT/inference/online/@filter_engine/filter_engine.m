function engine = filter_engine(tbn_engine)
% FILTER_ENGINE Create an engine which does online filtering
% function engine = filter_engine(tbn_engine)

engine.tbn_engine = tbn_engine;
engine.f = [];  % space to store filtered message
engine.old_f = [];
engine.b = []; % space to store smoothed message
engine.t = [];
engine = class(engine, 'filter_engine');
