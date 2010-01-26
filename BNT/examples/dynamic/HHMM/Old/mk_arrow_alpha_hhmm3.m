% Make the following HHMM
%
%     LH                  RH
%    /                      \
%   /                        \
%  LR -> UD -> RL -> DU       RL -> UD -> LR -> DU
%   \
%    \
%     Q1 -> Q2
%
% where level 1 is fully interconnected (not shown)
% level 2 is left-right
% and each model at level 3 is a 2 state LR shared HMM 

Qsizes = [2 4 2];
D = 3;

% LEVEL 1

startprob1 = 'ergodic';
transprob1 = 'ergodic';


% LEVEL 2

startprob = zeros(2, 4);
%        Q1  Q2
startprob(1, 1) = 1;
startprob(2, 3) = 1;

transprob = zeros(2, 4, 4);
transprob(1,:,:) = [0 1 0 0
		    0 0 1 0
		    0 0 0 1
		    0 0 0 1];
transprob(2,:,:) = [0 0 0 1
		    1 0 0 0
		    0 1 0 0
		    0 0 0 1];

Q2args = {'startprob', startprob, 'transprob', transprob};

% always terminate in state 4 (default)
% F2args

% LEVEL 3

% Defaults are fine: always start in state 1, left-right model, finish in state 2


% OBS LEVEl

chars = ['L', 'l', 'U', 'u', 'R', 'r', 'D', 'd'];
Osize = length(chars);

obsprob = zeros([4 2 Osize]);
%       Q2 Q3 O
obsprob(1, 1, find(chars == 'L')) =  1.0;
obsprob(1, 2, find(chars == 'l')) =  1.0;

obsprob(2, 1, find(chars == 'U')) =  1.0;
obsprob(2, 2, find(chars == 'u')) =  1.0;

obsprob(3, 1, find(chars == 'R')) =  1.0;
obsprob(3, 2, find(chars == 'r')) =  1.0;

obsprob(4, 1, find(chars == 'D')) =  1.0;
obsprob(4, 2, find(chars == 'd')) =  1.0;

Oargs = {'CPT', obsprob};


bnet = mk_hhmm3('Qsizes', Qsizes, 'Osize', Osize', 'discrete_obs', 1, 'Oargs', Oargs, 'Q1args', Q1args, 'Q2args', Q2args);

T = 20;
usecell = 0;
evidence = sample_dbn(bnet, T, usecell);      
%chars(evidence(end,:))

Q1 = 1; Q2 = 2; Q3 = 3; F3 = 4; F2 = 5; obs = 6;
Qnodes = [Q1 Q2 Q3]; Fnodes = [F2 F3];

pretty_print_hhmm_parse(evidence, Qnodes, Fnodes, obs, chars);

eclass = bnet.equiv_class;
S=struct(bnet.CPD{eclass(Q2,2)})
