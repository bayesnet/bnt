function bnet = mk_hhmm3(varargin)
% MK_HHMM3 Make a 3 level Hierarchical HMM
% bnet = mk_hhmm3(...)
%
% 3-layer hierarchical HMM where level 1 only connects to level 2, not 3 or obs.
% This enforces sub-models (which differ only in their Q1 index) to be shared.
% Also, we enforce the fact that each model always starts in its initial state
% and only finishes in its final state. However, the prob. of finishing (as opposed to
% self-transitioning to the final state) can be learned.
% The fact that we always finish from the same state means we do not need to condition
% F(i) on Q(i-1), since finishing prob is indep of calling context.
%
% The DBN is the same as Fig 10 in my tech report.
%
%   Q1 ---------->  Q1
%   |              / |
%   |             /  |
%   |  F2 -------    |
%   |  ^         \   |
%   | /|          \  |
%   v  |           v v
%   Q2-| -------->   Q2
%  /|  |             ^
% / |  |            /|
% | |  F3 ---------/ |
% | |  ^           \ |
% | v /              v
% | Q3 ----------->  Q3
% |  |    
% \  | 
%  v v    
%   O
%
% Q1 (slice 1) is clamped to be uniform.
% Q2 (slice 1) is clamped to always start in state 1.
% Q3 (slice 1) is clamped to always start in state 1.
% F3 by default will only finish if Q3 is in its last state (F3 is a tabular_CPD)
% F2 by default gets the default hhmmF_CPD params.
% Q1:Q3 (slice 2) by default gets the default hhmmQ_CPD params.
% O by default gets the default tabular/Gaussian params.
%
% Optional arguments in name/value format [default]
%
% Qsizes      - sizes at each level [ none ]
% Osize       - size of O node [ none ]
% discrete_obs - 1 means O is tabular_CPD, 0 means O is gaussian_CPD [0]
% Oargs       - cell array of args to pass to the O CPD  [ {} ]
% Q1args      - args to be passed to constructor for Q1 (slice 2) [ {} ]
% Q2args      - args to be passed to constructor for Q2 (slice 2) [ {} ]
% Q3args      - args to be passed to constructor for Q3 (slice 2) [ {} ]
% F2args       - args to be passed to constructor for F2 [ {} ]
% F3args       - args to be passed to constructor for F3 [ {'CPT', finish in last Q3 state} ]
%

ss = 6; D = 3;
Q1 = 1; Q2 = 2; Q3 = 3; F3 = 4; F2 = 5; obs = 6;
Qnodes = [Q1 Q2 Q3]; Fnodes = [F2 F3];
names = {'Q1', 'Q2', 'Q3', 'F3', 'F2', 'obs'};

intra = zeros(ss);
intra(Q1, Q2) = 1;
intra(Q2, [F2 Q3 obs]) = 1;
intra(Q3, [F3 obs]) = 1;
intra(F3, F2) = 1;

inter = zeros(ss);
inter(Q1,Q1) = 1;
inter(Q2,Q2) = 1;
inter(Q3,Q3) = 1;
inter(F2,[Q1 Q2]) = 1;
inter(F3,[Q2 Q3]) = 1;


% get sizes of nodes
args = varargin;
nargs = length(args);
Qsizes = [];
Osize = 0;
for i=1:2:nargs
  switch args{i},
   case 'Qsizes', Qsizes = args{i+1}; 
   case 'Osize', Osize = args{i+1}; 
  end
end
if isempty(Qsizes), error('must specify Qsizes'); end
if Osize==0, error('must specify Osize'); end
  
% set default params
discrete_obs = 0;
Oargs = {};
Q1args = {};
Q2args = {};
Q3args = {};
F2args = {};

% P(Q3, F3)
CPT = zeros(Qsizes(3), 2);
% Each model can only terminate in its final state.
% 0 params will remain 0 during EM, thus enforcing this constraint.
CPT(:, 1) = 1.0; % all states turn F off ...
p = 0.5;
CPT(Qsizes(3), 2) = p; % except the last one
CPT(Qsizes(3), 1) = 1-p;
F3args = {'CPT', CPT};

for i=1:2:nargs
  switch args{i},
   case 'discrete_obs', discrete_obs = args{i+1}; 
   case 'Oargs',        Oargs = args{i+1};
   case 'Q1args',       Q1args = args{i+1};
   case 'Q2args',       Q2args = args{i+1};
   case 'Q3args',       Q3args = args{i+1};
   case 'F2args',       F2args = args{i+1};
   case 'F3args',       F3args = args{i+1};
  end
end

ns = zeros(1,ss);
ns(Qnodes) = Qsizes;
ns(obs) = Osize;
ns(Fnodes) = 2;

dnodes = [Qnodes Fnodes];
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
CPT = normalise(ones(1,ns(Q1)));
bnet.CPD{eclass(Q1,1)} = tabular_CPD(bnet, Q1, 'CPT', CPT, 'adjustable', 0);

% each model always starts in state 1
CPT = zeros(ns(Q1), ns(Q2));
CPT(:, 1) = 1.0;
bnet.CPD{eclass(Q2,1)} = tabular_CPD(bnet, Q2, 'CPT', CPT, 'adjustable', 0);

% each model always starts in state 1
CPT = zeros(ns(Q2), ns(Q3));
CPT(:, 1) = 1.0;
bnet.CPD{eclass(Q3,1)} = tabular_CPD(bnet, Q3, 'CPT', CPT, 'adjustable', 0);

bnet.CPD{eclass(F2,1)}  = hhmmF_CPD(bnet, F2, Qnodes, 2, D, F2args{:});

bnet.CPD{eclass(F3,1)}  = tabular_CPD(bnet, F3, F3args{:});

if discrete_obs
  bnet.CPD{eclass(obs,1)} = tabular_CPD(bnet, obs, Oargs{:});
else
  bnet.CPD{eclass(obs,1)} = gaussian_CPD(bnet, obs, Oargs{:});
end

% SLICE 2

bnet.CPD{eclass(Q1,2)} = hhmmQ_CPD(bnet, Q1+ss, Qnodes, 1, D, Q1args{:});
bnet.CPD{eclass(Q2,2)} = hhmmQ_CPD(bnet, Q2+ss, Qnodes, 2, D, Q2args{:});
bnet.CPD{eclass(Q3,2)} = hhmmQ_CPD(bnet, Q3+ss, Qnodes, 3, D, Q3args{:});
