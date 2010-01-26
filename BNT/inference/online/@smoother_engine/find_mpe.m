function mpe = find_mpe(engine, ev)
% FIND_MPE Find the most probable explanation (Viterbi)
% mpe = enter_evidence(engine, evidence, ...)
%
% evidence{i,t} = [] if if X(i,t) is hidden, and otherwise contains its observed value (scalar or column vector)
%

mpe = cell(size(ev));
engine.tbn_engine = set_fields(engine.tbn_engine, 'maximize', 1);

T = size(ev, 2);
f = cell(1,T);
b = cell(1,T); % b{t}.clpot{c}
ll = zeros(1,T);
[f{1}, ll(1)] = fwd1(engine.tbn_engine, ev(:,1), 1);
for t=2:T
  [f{t}, ll(t)] = fwd(engine.tbn_engine, f{t-1}, ev(:,t), t);
end

if T==1
  [b{1}, mpe(:,1)] = backT_mpe(engine.tbn_engine, f{1}, ev(:,1), 1);
else
  [b{T}, mpe(:,T)] = backT_mpe(engine.tbn_engine, f{T}, ev(:,T-1:T), T);
  for t=T-1:-1:2
    [b{t}, mpe(:,t)] = back_mpe(engine.tbn_engine, b{t+1}, f{t}, ev(:,t-1:t), t);
  end
  t = 1;
  [b{t}, mpe(:,t)] = back1_mpe(engine.tbn_engine, b{t+1}, f{t}, ev(:,1), t);
end
engine.b = b;
  


