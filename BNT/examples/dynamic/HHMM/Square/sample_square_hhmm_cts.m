% Generate samples from the HHMM with the true params.

seed = 1;
rand('state', seed);
randn('state', seed);

discrete_obs = 0;

bnet = mk_square_hhmm(discrete_obs, 1);
Q1 = 1; Q2 = 2; Q3 = 3; F3 = 4; F2 = 5; Onode = 6;
Qnodes = [Q1 Q2 Q3]; Fnodes = [F2 F3];

for seqi=1:1
  evidence = sample_dbn(bnet, 'stop_test', 'is_F2_true_D3');      
  clf
  plot_square_hhmm(evidence);
  %pretty_print_hhmm_parse(evidence, Qnodes, Fnodes, Onode, []);
  fprintf('sequence %d has length %d; press key to continue\n', seqi, size(evidence,2))
  pause
end
