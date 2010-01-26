% Compare online filtering algorithms on some DBNs

seed = 0;
rand('state', seed);
randn('state', seed);

if 0
  N = 3;
  Q = 2;
  obs_size = 1;
  discrete_obs = 0;
  bnet = mk_chmm(N, Q, obs_size, discrete_obs);
else
  %bnet = mk_bat_dbn;
  bnet = mk_water_dbn;
end

T = 3;

engine = {};
engine{end+1} = filter_engine(hmm_2TBN_inf_engine(bnet));
engine{end+1} = filter_engine(jtree_2TBN_inf_engine(bnet));

time = cmp_online_inference(bnet, engine, T);
