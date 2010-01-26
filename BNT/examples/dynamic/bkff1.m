% Compare different implementations of fully factored Boyen Koller

water = 1;
if water
  bnet = mk_water_dbn;
else
  N = 5;
  Q = 2;
  Y = 2;
  bnet = mk_chmm(N, Q, Y);
end
ss = length(bnet.intra);

engine = {};
engine{end+1} = bk_inf_engine(bnet, 'clusters', 'ff');        
engine{end+1} = bk_ff_hmm_inf_engine(bnet);   
E = length(engine);

T = 5;
time = cmp_inference_dbn(bnet, engine, T, 'singletons_only', 1)



