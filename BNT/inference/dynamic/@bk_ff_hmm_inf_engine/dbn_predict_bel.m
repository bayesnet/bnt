function engine = dbn_predict_bel(engine, lag)
% DBN_PREDICT_BEL Predict the belief state 'lag' steps into the future (bk_ff_hmm)
% engine = dbn_predict_bel(engine, lag)
% 'lag' defaults to 1

if nargin < 2, lag = 1; end

for d=1:lag
  %newbel = engine.transmat' * engine.bel;
  newbel = normalise(engine.transmat' * engine.bel); 
  
  hnodes = engine.hnodes;
  bnet = bnet_from_engine(engine);
  ns = bnet.node_sizes;
  [marginals, marginalsT] = project_joint_onto_marginals(newbel, hnodes, ns);
  newbel = combine_marginals_into_joint(marginalsT, hnodes, ns);          
  engine.bel_marginals = marginalsT;
  engine.bel = newbel;
end
