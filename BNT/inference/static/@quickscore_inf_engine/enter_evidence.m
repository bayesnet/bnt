function engine = enter_evidence(engine, pos, neg)
% ENTER_EVIDENCE Add evidence to the QMR network
% engine = enter_evidence(engine, pos, neg)
%
% pos = list of leaves that have positive observations
% neg = list of leaves that have negative observations

% Extract params for the observed findings
obs = myunion(pos, neg);
%inhibit_obs = engine.inhibit(obs, :);
inhibit_obs = engine.inhibit(:,obs)';
leak_obs = engine.leak(obs);

% Find what nodes correspond to the original observed leaves
pos2 = find_equiv_posns(pos, obs);
neg2 = find_equiv_posns(neg, obs);
engine.post = quickscore(pos2, neg2, inhibit_obs, engine.prior, leak_obs); 
%engine.post = C_quickscore(pos2, neg2, inhibit_obs, engine.prior, leak_obs); 

