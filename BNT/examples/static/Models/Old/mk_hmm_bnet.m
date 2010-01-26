function [bnet, onodes] = mk_hmm_bnet(T, Q, O, cts_obs, param_tying)
% MK_HMM_BNET Make a (static( bnet to represent a hidden Markov model
% [bnet, onodes] = mk_hmm_bnet(T, Q, O, cts_obs, param_tying)
%
% T = num time slices
% Q = num hidden states
% O = size of the observed node (num discrete values or length of vector)
% cts_obs - 1 means the observed node is a continuous-valued vector, 0 means it's discrete
% param_tying - 1 means we create 3 CPDs, 0 means we create 1 CPD per node

N = 2*T;
dag = zeros(N);
for i=1:T-1
  dag(i,i+1)=1;
end
onodes = T+1:N;
for i=1:T
  dag(i, onodes(i)) = 1;
end

if cts_obs
  dnodes = 1:T;
else
  dnodes = 1:N;
end
ns = [Q*ones(1,T) O*ones(1,T)];

if param_tying
  eclass = [1 2*ones(1,T-1) 3*ones(1,T)];
else
  eclass = 1:N;
end

bnet = mk_bnet(dag, ns, dnodes, eclass);

hnodes = mysetdiff(1:N, onodes);
if ~param_tying
  for i=hnodes(:)'
    bnet.CPD{i} = tabular_CPD(bnet, i);
  end
  if cts_obs
    for i=onodes(:)'
      bnet.CPD{i} = gaussian_CPD(bnet, i);
    end
  else
    for i=onodes(:)'
      bnet.CPD{i} = tabular_CPD(bnet, i);
    end
  end
else
  bnet.CPD{1} = tabular_CPD(bnet, 1);
  bnet.CPD{2} = tabular_CPD(bnet, 2);
  if cts_obs
    bnet.CPD{3} = gaussian_CPD(bnet, 3);
  else
    bnet.CPD{3} = tabular_CPD(bnet, 3);
  end
end
