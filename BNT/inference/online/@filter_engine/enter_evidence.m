function [engine, LL] = enter_evidence(engine, ev, t)
% ENTER_EVIDENCE Call the online filter
% [engine, loglik] = enter_evidence(engine, evidence, ...)
%
% evidence{i} = [] if if X(i) is hidden, and otherwise contains its observed value (scalar or column vector)

engine.old_f = engine.f;
if t==1
  [engine.f, LL] = fwd1(engine.tbn_engine, ev, 1);
else
  [engine.f, LL] = fwd(engine.tbn_engine, engine.old_f, ev, t);
end
engine.b = backT(engine.tbn_engine, engine.f, t);
engine.t = t;
