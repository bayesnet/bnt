function [bnet, Qnodes, Fnodes, Onode] = mk_hhmm(varargin)
% MK_HHMM Make a Hierarchical HMM
% function [bnet, Qnodes, Fnodes, Onode] = mk_hhmm(...)
%
% e.g. 3-layer hierarchical HMM where level 1 only connects to level 2
% and the parents of the observed node are levels 2 and 3.
% (This DBN is the same as Fig 10 in my tech report.)
%
%   Q1 ---------->   Q1
%   |  \            ^ |
%   |   v          /  |
%   |    F2 ------/   |
%   |   ^ ^        \  |
%   |  /  |         \ |
%   | /   |          ||
%   v     |          vv
%   Q2----| --------> Q2
%  /| \   |          ^|
% / |  v  |         / |
% | |   F3 --------/  |
% | |   ^          \  |
% | v  /            v v
% | Q3 ----------->  Q3
% |  |    
% \  | 
%  v v    
%   O
%
%
% Optional arguments in name/value format [default value in brackets]
%
% Qsizes      - sizes at each level [ none ]
% allQ       - 1 means level i connects to all Q levels below, 0 means just to i+1 [0]
% transprob  - transprob{d}(i,k,j) = P(Q(d,t)=j|Q(d,t-1)=i,Q(1:d-1,t)=k)  ['leftright']
% startprob  - startprob{d}(k,j) = P(Q(d,t)=j|Q(1:d-1,t)=k)  ['leftstart']
% termprob   - termprob{d}(k,j) = P(F(d,t)=2|Q(1:d-1,t)=k,Q(d,t)=j) for d>1 ['rightstop']
% selfprop   - prob of a self transition (termprob default = 1-selfprop) [0.8]
% Osize       - size of O node
% discrete_obs - 1 means O is tabular_CPD, 0 means gaussian_CPD [0]
% Oargs       - cell array of args to pass to the O CPD  [ {} ]
% Ops         - Q parents of O [Qnodes(end)]
% F1          - 1 means level 1 can finish (restart), else there is no F1->Q1 arc [0]
% clamp1     - 1 means we clamp the params of the Q nodes in slice 1 (Qt1params) [1]
%   Note: the Qt1params are startprob, which should be shared with other slices.
%   However, in the current implementation, the Qt1params will only be estimated
%   from the initial state of each sequence.
%
% For d=1, startprob{1}(1,j) is only used in the first slice and
% termprob{1} is ignored, since we assume the top level never resets.
% Also, transprob{1}(i,j) can be used instead of transprob{1}(i,1,j).
%
% leftstart means the model always starts in state 1.
% rightstop means the model can only finish in its last state (Qsize(d)).
% unif means each state is equally like to reach any other
% rnd means the transition/starting probs are random (drawn from rand)
%
% Q1:QD in slice 1 are of type tabular_CPD
% Q1:QD in slice 2 are of type hhmmQ_CPD.
% F(2:D-1) is of type hhmmF_CPD, FD is of type tabular_CPD.

args = varargin;
nargs = length(args);

% get sizes of nodes and topology
Qsizes = [];
Osize = [];
allQ = 0;
Ops = [];
F1 = 0;
for i=1:2:nargs
  switch args{i},
   case 'Qsizes', Qsizes = args{i+1}; 
   case 'Osize',  Osize = args{i+1}; 
   case 'allQ',   allQ = args{i+1}; 
   case 'Ops',    Ops = args{i+1}; 
   case 'F1',     F1 = args{i+1}; 
  end
end
if isempty(Qsizes), error('must specify Qsizes'); end
if Osize==0, error('must specify Osize'); end
D = length(Qsizes);
Qnodes = 1:D;

if isempty(Ops), Ops = Qnodes(end); end


[intra, inter, Qnodes, Fnodes, Onode] = mk_hhmm_topo(D, allQ, Ops, F1);
ss = length(intra);
names = {};

if F1
  Fnodes_ndx = Fnodes;
else
  Fnodes_ndx = [-1 Fnodes]; % Fnodes(1) is a dummy index
end
		
% set default params
discrete_obs = 0;
Oargs = {};
startprob = cell(1,D);
startprob{1} = 'unif';
for d=2:D
  startprob{d} = 'leftstart';
end
transprob = cell(1,D);
transprob{1} = 'unif';
for d=2:D
  transprob{d} = 'leftright';
end
termprob = cell(1,D);
for d=2:D
  termprob{d} = 'rightstop';
end
selfprob = 0.8;
clamp1 = 1;

for i=1:2:nargs
  switch args{i},
   case 'discrete_obs', discrete_obs = args{i+1}; 
   case 'Oargs',        Oargs = args{i+1};
   case 'startprob',    startprob = args{i+1};
   case 'transprob',    transprob = args{i+1};
   case 'termprob',     termprob = args{i+1};
   case 'selfprob',     selfprob = args{i+1};
   case 'clamp1',       clamp1 = args{i+1};
  end
end

ns = zeros(1,ss);
ns(Qnodes) = Qsizes;
ns(Onode) = Osize;
ns(Fnodes) = 2;

dnodes = [Qnodes Fnodes];
if discrete_obs
  dnodes = [dnodes Onode];
end
onodes = [Onode];

bnet = mk_dbn(intra, inter, ns, 'observed', onodes, 'discrete', dnodes, 'names', names);
eclass = bnet.equiv_class;

for d=1:D
  if d==1
    Qps = [];
  elseif allQ
    Qps = Qnodes(1:d-1);
  else
    Qps = Qnodes(d-1);
  end
  Qpsz = prod(ns(Qps));
  Qsz = ns(Qnodes(d));
  if isstr(startprob{d})
    switch startprob{d}
     case 'unif', startprob{d} = mk_stochastic(ones(Qpsz, Qsz));
     case 'rnd', startprob{d} = mk_stochastic(rand(Qpsz, Qsz));
     case 'leftstart', startprob{d} = zeros(Qpsz, Qsz); startprob{d}(:,1) = 1;
    end
  end
  if isstr(transprob{d})
    switch transprob{d}
     case 'unif', transprob{d} = mk_stochastic(ones(Qsz, Qpsz, Qsz));
     case 'rnd', transprob{d} = mk_stochastic(rand(Qsz, Qpsz, Qsz));
     case 'leftright',
      LR = mk_leftright_transmat(Qsz, selfprob);
      temp = repmat(reshape(LR, [1 Qsz Qsz]), [Qpsz 1 1]); % transprob(k,i,j)
      transprob{d} = permute(temp, [2 1 3]); % now transprob(i,k,j)
    end
  end
  if isstr(termprob{d})
    switch termprob{d}
     case 'unif', termprob{d} = mk_stochastic(ones(Qpsz, Qsz, 2));
     case 'rnd', termprob{d} = mk_stochastic(rand(Qpsz, Qsz, 2));
     case 'rightstop',
      %termprob(k,i,t) Might terminate if i=Qsz; will not terminate if i<Qsz
      stopprob = 1-selfprob;
      termprob{d} = zeros(Qpsz, Qsz, 2);
      termprob{d}(:,Qsz,2) = stopprob;
      termprob{d}(:,Qsz,1) = 1-stopprob;
      termprob{d}(:,1:(Qsz-1),1) = 1;
     otherwise, error(['unrecognized termprob ' termprob{d}])
    end
  elseif d>1 % passed in termprob{d}(k,j)
    temp = termprob{d};
    termprob{d} = zeros(Qpsz, Qsz, 2);
    termprob{d}(:,:,2) = temp;
    termprob{d}(:,:,1) = ones(Qpsz,Qsz) - temp;
  end
end


% SLICE 1

for d=1:D
  bnet.CPD{eclass(Qnodes(d),1)} = tabular_CPD(bnet, Qnodes(d), 'CPT', startprob{d}, 'adjustable', clamp1);
end

if F1
  d = 1;
  bnet.CPD{eclass(Fnodes_ndx(d),1)}  = hhmmF_CPD(bnet, Fnodes_ndx(d), Qnodes(d), Fnodes_ndx(d+1), ...
						 'termprob', termprob{d});
end
for d=2:D-1
  if allQ
    Qps = Qnodes(1:d-1);
  else
    Qps = Qnodes(d-1);
  end
  bnet.CPD{eclass(Fnodes_ndx(d),1)}  = hhmmF_CPD(bnet, Fnodes_ndx(d), Qnodes(d), Fnodes_ndx(d+1), ...
						 'Qps', Qps, 'termprob', termprob{d});
end
bnet.CPD{eclass(Fnodes_ndx(D),1)}  = tabular_CPD(bnet, Fnodes_ndx(D), 'CPT', termprob{D});

if discrete_obs
  bnet.CPD{eclass(Onode,1)} = tabular_CPD(bnet, Onode, Oargs{:});
else
  bnet.CPD{eclass(Onode,1)} = gaussian_CPD(bnet, Onode, Oargs{:});
end

% SLICE 2

%for d=1:D
%  bnet.CPD{eclass(Qnodes(d),2)} = hhmmQ_CPD(bnet, Qnodes(d)+ss, Qnodes, d, D, ...
%					    'startprob', startprob{d}, 'transprob', transprob{d}, ...
%					    'allQ', allQ);
%end

d = 1;
if F1
  bnet.CPD{eclass(Qnodes(d),2)} = hhmmQ_CPD(bnet, Qnodes(d)+ss, 'Fself', Fnodes_ndx(d), ...
					    'Fbelow', Fnodes_ndx(d+1), ...
					    'startprob', startprob{d}, 'transprob', transprob{d});
else
    bnet.CPD{eclass(Qnodes(d),2)} = hhmmQ_CPD(bnet, Qnodes(d)+ss, ...
					    'Fbelow', Fnodes_ndx(d+1), ...
					    'startprob', startprob{d}, 'transprob', transprob{d});
end
for d=2:D-1
  if allQ
    Qps = Qnodes(1:d-1);
  else
    Qps = Qnodes(d-1);
  end
  Qps = Qps + ss; % since all in slice 2
  bnet.CPD{eclass(Qnodes(d),2)} = hhmmQ_CPD(bnet, Qnodes(d)+ss, 'Fself', Fnodes_ndx(d), ...
					    'Fbelow', Fnodes_ndx(d+1), 'Qps', Qps, ...
					    'startprob', startprob{d}, 'transprob', transprob{d});
end
d = D;
if allQ
  Qps = Qnodes(1:d-1);
else
  Qps = Qnodes(d-1);
end
Qps = Qps + ss; % since all in slice 2
bnet.CPD{eclass(Qnodes(d),2)} = hhmmQ_CPD(bnet, Qnodes(d)+ss, 'Fself', Fnodes_ndx(d), ...
					  'Qps', Qps, ...
					  'startprob', startprob{d}, 'transprob', transprob{d});
