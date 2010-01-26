% Make a DBN with the following inter-connectivity matrix
%    1
%   /  \
%  2   3
%   \ /
%    4 
%    |
%    5
% where all arcs point down. In addition, there are persistence arcs from each node to itself.
% There are no intra-slice connections.
% Nodes have noisy-or CPDs.
% Node 1 turns on spontaneously due to its leaky source.
% This effect trickles down to the other nodes in the order shown.
% All the other nodes inhibit their leaks.
% None of the nodes inhibit the connection from themselves, so that once they are on, they remain
% on (persistence).
%
% This model was used in the experiments reported in
% - "Learning the structure of DBNs", Friedman, Murphy and Russell, UAI 1998.
% where the structure was learned even in the presence of missing data.
% In that paper, we used the structural EM algorithm.
% Here, we assume full observability and tabular CPDs for the learner, so we can use a much
% simpler learning algorithm.

ss = 5;

inter = eye(ss);
inter(1,[2 3]) = 1;
inter(2,4)=1;
inter(3,4)=1;
inter(4,5)=1;

intra = zeros(ss);
ns = 2*ones(1,ss);

bnet = mk_dbn(intra, inter, ns);

% All nodes start out off
for i=1:ss
  bnet.CPD{i} = tabular_CPD(bnet, i, [1.0 0.0]');
end

% The following params correspond to Fig 4a in the UAI 98 paper
% The first arg is the leak inhibition prob.
% The vector contains the inhib probs from the parents in the previous slice;
% the last element is self, which is never inhibited.
bnet.CPD{1+ss} = noisyor_CPD(bnet, 1+ss, 0.8, 0);
bnet.CPD{2+ss} = noisyor_CPD(bnet, 2+ss, 1, [0.9 0]);
bnet.CPD{3+ss} = noisyor_CPD(bnet, 3+ss, 1, [0.8 0]);
bnet.CPD{4+ss} = noisyor_CPD(bnet, 4+ss, 1, [0.7 0.6 0]);
bnet.CPD{5+ss} = noisyor_CPD(bnet, 5+ss, 1, [0.5 0]);


% Generate some training data

nseqs = 20;
seqs = cell(1,nseqs);
T = 30;
for i=1:nseqs
  seqs{i} = sample_dbn(bnet, T);
end

max_fan_in = 3; % let's cheat a little here

% computing num. incorrect edges as a fn of the size of the training set
%sz = [5 10 15 20];    
sz = [5 10];    
h = zeros(1, length(sz));
for i=1:length(sz)
  inter2 = learn_struct_dbn_reveal(seqs(1:sz(i)), ns, max_fan_in);
  h(i) = sum(abs(inter(:)-inter2(:))); % hamming distance
end
h
