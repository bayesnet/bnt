function engine = smoother_engine(tbn_engine)
% SMOOTHER_ENGINE Create an engine which does offline (fixed-interval) smoothing in O(T) space/time
% function engine = smoother_engine(tbn_engine)
%
% tbn_engine is any 2TBN inference engine which supports the following methods:
% fwd, fwd1, back, backT, back, marginal_nodes and marginal_family.

engine.tbn_engine = tbn_engine;
engine.b = []; % space to store smoothed messages
engine = class(engine, 'smoother_engine');
%engine = class(engine, 'smoother_engine', inf_engine(bnet_from_engine(tbn_engine)));

