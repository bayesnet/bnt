function bnet = mk_rnd_map_hhmm(varargin)

% We copy the deterministic structure of the real HHMM,
% but randomize the probabilities of the adjustable CPDs.
% The key trick is that 0s in the real HHMM remain 0
% even when multiplied by a randon number.

obs_model = 'unique';

for i=1:2:length(varargin)
  switch varargin{i},
   case 'obs_model', obs_model = varargin{i+1};
  end
end


unique_obs = strcmp(obs_model, 'unique');

psuccess = 0.9;
% must be less than 1, so that pfail > 0
% otherwise we copy too many 0s
bnet = mk_map_hhmm('p', psuccess, 'obs_model', obs_model);
ns = bnet.node_sizes;
ss = bnet.nnodes_per_slice;

U = 1; A = 2; C = 3; F = 4;
%unique_obs = (bnet.nnodes_per_slice == 5);
if unique_obs
  onodes = 5;
else
  north = 5; east = 6; south = 7; west = 8;
  onodes = [north east south west];
end

eclass = bnet.equiv_class;
S=struct(bnet.CPD{eclass(F,1)});
CPT = mk_stochastic(rand(size(S.CPT)) .* S.CPT);
bnet.CPD{eclass(F,1)} = tabular_CPD(bnet, F, 'CPT', CPT);


% Observation model
if unique_obs
  CPT = zeros(ns(A), ns(C), 5);
  CPT(1,1,1)=1;  % Theo state 4
  CPT(1,2,2)=1;  % Theo state 5
  CPT(1,3,3)=1; % Theo state 6
  CPT(2,1,4)=1; % Theo state 9
  CPT(2,2,5)=1; % Theo state 10
  %CPT(2,3,:) undefined
  O = onodes(1);
  bnet.CPD{eclass(O,1)} = tabular_CPD(bnet, O, 'CPT', CPT);
else
  for i=[north east south west]
    CPT = mk_stochastic(rand(ns(A), ns(C), 2));
    bnet.CPD{eclass(i,1)} = tabular_CPD(bnet, i, 'CPT', CPT);
  end
end

% Define the CPDs for slice 2

startprob = zeros(ns(U), ns(A));
S = struct(bnet.CPD{eclass(A,2)});
transprob = mk_stochastic(rand(size(S.transprob)) .* S.transprob);
bnet.CPD{eclass(A,2)} = hhmm2Q_CPD(bnet, A+ss, 'Fbelow', F, ...
				  'startprob', startprob, 'transprob', transprob);

S = struct(bnet.CPD{eclass(C,2)});
transprob = mk_stochastic(rand(size(S.transprob)) .* S.transprob);
startprob = mk_stochastic(rand(size(S.startprob)) .* S.startprob);
bnet.CPD{eclass(C,2)} = hhmm2Q_CPD(bnet, C+ss, 'Fself', F, ...
				  'startprob', startprob, 'transprob', transprob);


