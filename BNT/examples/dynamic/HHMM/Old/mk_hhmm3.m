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
%
% Optional arguments in name/value format [default]
%
% Qsizes      - sizes at each level [ none ]
% Osize       - size of O node [ none ]
% discrete_obs - 1 means O is tabular_CPD, 0 means O is gaussian_CPD [0]
% Oargs       - cell array of args to pass to the O CPD  [ {} ]
% transprob1  - transprob1(i,j) = P(Q1(t)=j|Q1(t-1)=i)  ['ergodic']
% startprob1  - startprob1(j) = P(Q1(t)=j)  ['leftstart']
% transprob2  - transprob2(i,k,j) = P(Q2(t)=j|Q2(t-1)=i,Q1(t)=k)  ['leftright']
% startprob2  - startprob2(k,j) = P(Q2(t)=j|Q1(t)=k)  ['leftstart']
% termprob2   - termprob2(j,f) = P(F2(t)=f|Q2(t)=j)  ['rightstop']
% transprob3  - transprob3(i,k,j) = P(Q3(t)=j|Q3(t-1)=i,Q2(t)=k)  ['leftright']
% startprob3  - startprob3(k,j) = P(Q3(t)=j|Q2(t)=k)  ['leftstart']
% termprob3   - termprob3(j,f) = P(F3(t)=f|Q3(t)=j)  ['rightstop']
%
% leftstart means the model always starts in state 1.
% rightstop means the model always finished in its last state (Qsize(d)).
%
% Q1:Q3 in slice 1 are of type tabular_CPD
% Q1:Q3 in slice 2 are of type hhmmQ_CPD.
% F2 is of type hhmmF_CPD, F3 is of type tabular_CPD.

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
startprob1 = 'ergodic';
startprob2 = 'leftstart';
startprob3 = 'leftstart';
transprob1 = 'ergodic';
transprob2 = 'leftright';
transprob3 = 'leftright';
termprob2 = 'rightstop';
termprob3 = 'rightstop';


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

if strcmp(startprob1, 'ergodic')
  startprob1 = normalise(ones(1,ns(Q1)));
end
if strcmp(startprob2, 'leftstart')
  startprob2 = zeros(ns(Q1), ns(Q2));
  starpbrob2(:, 1) = 1.0;
end
if strcmp(startprob3, 'leftstart')
  startprob3 = zeros(ns(Q2), ns(Q3));
  starpbrob3(:, 1) = 1.0;
end

if strcmp(termprob2, 'rightstop')
  p = 0.9;
  termprob2 = zeros(Qsize(2),2);
  termprob2(:, 2) = p; 
  termprob2(:, 1) = 1-p; 
  termprob2(1:(Qsize(2)-1), 1) = 1; 
end
if strcmp(termprob3, 'rightstop')
  p = 0.9;
  termprob3 = zeros(Qsize(3),2);
  termprob3(:, 2) = p; 
  termprob3(:, 1) = 1-p; 
  termprob3(1:(Qsize(3)-1), 1) = 1; 
end


% SLICE 1

% We clamp untied nodes in the first slice, since their params can't be estimated
% from just one sequence

bnet.CPD{eclass(Q1,1)} = tabular_CPD(bnet, Q1, 'CPT', startprob1, 'adjustable', 0);
bnet.CPD{eclass(Q2,1)} = tabular_CPD(bnet, Q2, 'CPT', startprob2, 'adjustable', 0);
bnet.CPD{eclass(Q3,1)} = tabular_CPD(bnet, Q3, 'CPT', startprob3, 'adjustable', 0);

bnet.CPD{eclass(F2,1)}  = hhmmF_CPD(bnet, F2, Qnodes, 2, D, 'termprob', termprob2);
bnet.CPD{eclass(F3,1)}  = tabular_CPD(bnet, F3, 'CPT', termprob3);

if discrete_obs
  bnet.CPD{eclass(obs,1)} = tabular_CPD(bnet, obs, Oargs{:});
else
  bnet.CPD{eclass(obs,1)} = gaussian_CPD(bnet, obs, Oargs{:});
end

% SLICE 2

bnet.CPD{eclass(Q1,2)} = hhmmQ_CPD(bnet, Q1+ss, Qnodes, 1, D, 'transprob', transprob1, 'startprob', startprob1);
bnet.CPD{eclass(Q2,2)} = hhmmQ_CPD(bnet, Q2+ss, Qnodes, 2, D, 'transprob', transprob2, 'startprob', startprob2);
bnet.CPD{eclass(Q3,2)} = hhmmQ_CPD(bnet, Q3+ss, Qnodes, 3, D, 'transprob', transprob3, 'startprob', startprob3);

