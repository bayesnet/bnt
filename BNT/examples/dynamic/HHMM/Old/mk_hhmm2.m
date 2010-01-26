function bnet = mk_hhmm2(varargin)
% MK_HHMM2 Make a 2 level Hierarchical HMM
% bnet = mk_hhmm2(...)
%
% 2-layer hierarchical HMM  (node numbers in parens)
%
%   Q1(1) ---------> Q1(5)
% /  | \            / |
% |  |  v          /  |
% |  |  F2(3) --- /   |
% |  |  ^         \   |
% |  | /           \  |
% |  v              \ v
% |  Q2(2)--------> Q2 (6)
% |  |    
% \  | 
%  v v    
%   O(4)
%
%
% Optional arguments [default]
%
% discrete_obs - 1 means O is tabular_CPD, 0 means O is gaussian_CPD [0]
% obsCPT       - CPT(o,q1,q2) params for O ['rnd']
% mu           - mu(:,q1,q2) params for O [ [] ]
% Sigma        - Sigma(:,q1,q2) params for O [ [] ]
%
% F2toQ1       - 1 if Q2 is an hhmm_CPD, 0 if F2 -> Q2 arc is absent, so level 2 never resets [1]
% Q1args        - arguments to be passed to the constructors for Q1(t=2) [ {} ]
% Q2args        - arguments to be passed to the constructors for Q2(t=2) [ {} ]
%
% F2 only turns on (wp 0.5) when Q2 enters its final state.
% Q1 (slice 1) is clamped to be uniform.
% Q2 (slice 1) is clamped to always start in state 1.

[os nmodels nstates] = size(mu);

ss = 4;
Q1 = 1; Q2 = 2; F2 = 3; obs = 4;
Qnodes = [Q1 Q2];
names = {'Q1', 'Q2', 'F2', 'obs'};
intra = zeros(ss);
intra(Q1, [Q2 F2 obs]) = 1;
intra(Q2, [F2 obs]) = 1;

inter = zeros(ss);
inter(Q1,Q1) = 1;
inter(F2,Q1) = 1;
if F2toQ2
  inter(F2,Q2)=1;
end
inter(Q2,Q2) = 1;

ns = zeros(1,ss);

ns(Q1) = nmodels;
ns(Q2) = nstates;
ns(F2) = 2;
ns(obs) = os;

dnodes = [Q1 Q2 F2];
if discrete_obs
  dnodes = [dnodes obs];
end
onodes = [obs];

bnet = mk_dbn(intra, inter, ns, 'observed', onodes, 'discrete', dnodes, 'names', names);
eclass = bnet.equiv_class;

% SLICE 1

% We clamp untied nodes in the first slice, since their params can't be estimated
% from just one sequence

% uniform prior on initial model
CPT = normalise(ones(1,nmodels));
bnet.CPD{eclass(Q1,1)} = tabular_CPD(bnet, Q1, 'CPT', CPT, 'adjustable', 0);

% each model always starts in state 1
CPT = zeros(ns(Q1), ns(Q2));
CPT(:, 1) = 1.0;
bnet.CPD{eclass(Q2,1)} = tabular_CPD(bnet, Q2, 'CPT', CPT, 'adjustable', 0);

% Termination probability
CPT = zeros(ns(Q1), ns(Q2), 2);
if 1
  % Each model can only terminate in its final state.
  % 0 params will remain 0 during EM, thus enforcing this constraint.
  CPT(:, :, 1) = 1.0; % all states turn F off ...
  p = 0.5;
  CPT(:, ns(Q2), 2) = p; % except the last one
  CPT(:, ns(Q2), 1) = 1-p;
end
bnet.CPD{eclass(F2,1)}  = tabular_CPD(bnet, F2, 'CPT', CPT);

if discrete_obs
  bnet.CPD{eclass(obs,1)} = tabular_CPD(bnet, obs, obs_args{:});
else
  bnet.CPD{eclass(obs,1)} = gaussian_CPD(bnet, obs, obs_args{:});
end

% SLICE 2


bnet.CPD{eclass(Q1,2)} = hhmm_CPD(bnet, Q1+ss, Qnodes, 1, D, 'args', Q1args);

if F2toQ2
  bnet.CPD{eclass(Q2,2)} = hhmmQD_CPD(bnet, Q2+ss, Qnodes, 2, D, Q2args{:});
else
  bnet.CPD{eclass(Q2,2)} = tabular_CPD(bnet, Q2+ss, Q2args{:});
end
