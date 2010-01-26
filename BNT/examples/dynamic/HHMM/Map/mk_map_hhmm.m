function bnet = mk_map_hhmm(varargin)

% p is the prob of a successful move (defines the reliability of motors)
p = 1;
obs_model = 'unique';

for i=1:2:length(varargin)
  switch varargin{i},
   case 'p', p = varargin{i+1};
   case 'obs_model', obs_model = varargin{i+1};
  end
end


q = 1-p;
unique_obs = strcmp(obs_model, 'unique');

% assign numbers to the nodes in topological order
U = 1; A = 2; C = 3; F = 4;
if unique_obs
  onodes = 5;
else
  N = 5; E = 6; S = 7; W = 8; % north, east, south, west
  onodes = [N E S W];
end

% create graph structure

ss = 4 + length(onodes); % slice size
intra = zeros(ss,ss);
intra(U,F)=1;
intra(A,[C F onodes])=1;
intra(C,[F onodes])=1;

inter = zeros(ss,ss);
inter(U,[A C])=1;
inter(A,[A C])=1;
inter(F,[A C])=1;
inter(C,C)=1;

% node sizes
ns = zeros(1,ss);
ns(U) = 2; % left/right
ns(A) = 2;
ns(C) = 3;
ns(F) = 2;
if unique_obs
  ns(onodes) = 5; % we will assign each state a unique symbol
else
  ns(onodes) = 2;
end
l = 1; r = 2; % left/right
L = 1; R = 2;

% Make the DBN
bnet = mk_dbn(intra, inter, ns, 'observed', onodes);
eclass = bnet.equiv_class;



% Define CPDs for slice 1
% We clamp all the CPDs that are not tied,
% since we cannot learn them from a single sequence.

% uniform probs over actions (the input could be chosen from a policy)
bnet.CPD{eclass(U,1)} = tabular_CPD(bnet, U, 'CPT', mk_stochastic(ones(ns(U),1)), ...
				    'adjustable', 0);

% uniform probs over starting abstract state
bnet.CPD{eclass(A,1)} = tabular_CPD(bnet, A, 'CPT', mk_stochastic(ones(ns(A),1)), ...
				    'adjustable', 0);

% Uniform probs over starting concrete state, modulo the fact
% that corridor 2 is only of length 2.
CPT = zeros(ns(A), ns(C)); % CPT(i,j) = P(C starts in j | A=i)
CPT(1, :) = [1/3 1/3 1/3];
CPT(2, :) = [1/2 1/2 0];
bnet.CPD{eclass(C,1)} = tabular_CPD(bnet, C, 'CPT', CPT, 'adjustable', 0);

% Termination probs
CPT = zeros(ns(U), ns(A), ns(C), ns(F));
CPT(r,1,1,:) = [1 0];
CPT(r,1,2,:) = [1 0];
CPT(r,1,3,:) = [q p];
CPT(r,2,1,:) = [1 0];
CPT(r,2,2,:) = [q p];
CPT(l,1,1,:) = [q p];
CPT(l,1,2,:) = [1 0];
CPT(l,1,3,:) = [1 0];
CPT(l,2,1,:) = [q p];
CPT(l,2,2,:) = [1 0];

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
  % north/east/south/west can see wall (1) or opening (2)
  CPT = zeros(ns(A), ns(C), 2);
  CPT(:,:,1) = q;
  CPT(:,:,2) = p;
  bnet.CPD{eclass(W,1)} = tabular_CPD(bnet, W, 'CPT', CPT);
  bnet.CPD{eclass(E,1)} = tabular_CPD(bnet, E, 'CPT', CPT);
  CPT = zeros(ns(A), ns(C), 2);
  CPT(:,:,1) = p;
  CPT(:,:,2) = q;
  bnet.CPD{eclass(S,1)} = tabular_CPD(bnet, S, 'CPT', CPT);
  bnet.CPD{eclass(N,1)} = tabular_CPD(bnet, N, 'CPT', CPT);
end

% Define the CPDs for slice 2

% Abstract

% Since the top level never resets, the starting distribution is irrelevant:
% A2 will be determined by sampling from transmat(A1,:).
% But the code requires we specify it anyway; we make it all 0s, a dummy value.
startprob = zeros(ns(U), ns(A));

transmat = zeros(ns(U), ns(A), ns(A));
transmat(R,1,:) = [q p];
transmat(R,2,:) = [0 1];
transmat(L,1,:) = [1 0];
transmat(L,2,:) = [p q];

% Qps are the parents we condition the parameters on, in this case just
% the past action.
bnet.CPD{eclass(A,2)} = hhmm2Q_CPD(bnet, A+ss, 'Fbelow', F, ...
				  'startprob', startprob, 'transprob', transmat);



% Concrete

transmat = zeros(ns(C), ns(U), ns(A), ns(C));
transmat(1,r,1,:) = [q p 0.0];
transmat(2,r,1,:) = [0.0 q p];
transmat(3,r,1,:) = [0.0 0.0 1.0];
transmat(1,r,2,:) = [q p 0.0];
transmat(2,r,2,:) = [0.0 1.0 0.0];
%
transmat(1,l,1,:) = [1.0 0.0 0.0];
transmat(2,l,1,:) = [p q 0.0];
transmat(3,l,1,:) = [0.0 p q];
transmat(1,l,2,:) = [1.0 0.0 0.0];
transmat(2,l,2,:) = [p q 0.0];

% Add a new dimension for A(t-1), by copying old vals,
% so the matrix is the same size as startprob


transmat = reshape(transmat, [ns(C) ns(U) ns(A) 1 ns(C)]);
transmat = repmat(transmat, [1 1 1 ns(A) 1]);

% startprob(C(t-1), U(t-1), A(t-1), A(t), C(t))
startprob = zeros(ns(C), ns(U), ns(A), ns(A), ns(C));
startprob(1,L,1,1,:) = [1.0 0.0 0.0];
startprob(3,R,1,2,:) = [1.0 0.0 0.0];
startprob(3,R,1,1,:) = [0.0 0.0 1.0];
% 
startprob(1,L,2,1,:) = [0.0 0.0 010];
startprob(2,L,2,1,:) = [1.0 0.0 0.0];
startprob(2,R,2,2,:) = [0.0 1.0 0.0];

% want transmat(U,A,C,At,Ct), ie. in topo order
transmat = permute(transmat, [2 3 1 4 5]);
startprob  = permute(startprob, [2 3 1 4 5]);
bnet.CPD{eclass(C,2)} = hhmm2Q_CPD(bnet, C+ss, 'Fself', F, ...
				  'startprob', startprob, 'transprob', transmat);


