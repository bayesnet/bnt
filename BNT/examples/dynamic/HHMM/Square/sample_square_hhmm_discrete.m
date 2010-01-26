% Generate samples from the HHMM with the true params.

seed = 0;
rand('state', seed);
randn('state', seed);

discrete_obs = 1;

bnet = mk_square_hhmm(discrete_obs, 1);

Tmax = 30;
Q1 = 1; Q2 = 2; Q3 = 3; F3 = 4; F2 = 5; Onode = 6;
Qnodes = [Q1 Q2 Q3]; Fnodes = [F2 F3];
chars = ['L', 'l', 'U', 'u', 'R', 'r', 'D', 'd'];
  
for seqi=1:3
  evidence = cell2num(sample_dbn(bnet, 'stop_test', 'is_F2_true_D3'));
  T = size(evidence, 2)
  pretty_print_hhmm_parse(evidence, Qnodes, Fnodes, Onode, chars);
end
