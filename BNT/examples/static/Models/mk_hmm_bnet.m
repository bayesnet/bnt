function bnet = mk_hmm_bnet(T, Q, O, cts_obs, param_tying)
% MK_HMM_BNET Make a (static) bnet to represent a hidden Markov model
% bnet = mk_hmm_bnet(T, Q, O, cts_obs, param_tying)
%
% T = num time slices
% Q = num hidden states
% O = size of the observed node (num discrete values or length of vector)
% cts_obs - 1 means the observed node is a continuous-valued vector, 0 means it's discrete
% param_tying - 1 means we create 3 CPDs, 0 means we create 1 CPD per node

N = 2*T;
dag = zeros(N);
%hnodes = 1:2:2*T;
hnodes = 1:T;
for i=1:T-1
  dag(hnodes(i), hnodes(i+1))=1;
end
%onodes = 2:2:2*T;
onodes = T+1:2*T;
for i=1:T
  dag(hnodes(i), onodes(i)) = 1;
end

if cts_obs
  dnodes = hnodes;
else
  dnodes = 1:N;
end
ns = ones(1,N);
ns(hnodes) = Q;
ns(onodes) = O;

if param_tying
  H1class = 1; Hclass = 2; Oclass = 3;
  eclass = ones(1,N);
  eclass(hnodes(2:end)) = Hclass;
  eclass(hnodes(1)) = H1class;
  eclass(onodes) = Oclass;
else
  eclass = 1:N;
end

bnet = mk_bnet(dag, ns, 'observed', onodes, 'discrete', dnodes, 'equiv_class', eclass);

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
  bnet.CPD{H1class} = tabular_CPD(bnet, hnodes(1)); % prior
  bnet.CPD{Hclass} = tabular_CPD(bnet, hnodes(2)); % transition matrix
  if cts_obs
    bnet.CPD{Oclass} = gaussian_CPD(bnet, onodes(1));
  else
    bnet.CPD{Oclass} = tabular_CPD(bnet, onodes(1));
  end
end
