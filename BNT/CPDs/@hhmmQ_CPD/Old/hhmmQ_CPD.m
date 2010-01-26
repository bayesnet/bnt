function CPD = hhmmQ_CPD(bnet, self, Qnodes, d, D, varargin)
% HHMMQ_CPD Make the CPD for a Q node at depth D of a D-level hierarchical HMM
% CPD = hhmmQ_CPD(bnet, self, Qnodes, d, D, ...)
%
%  Fd(t-1) \   Q1:d-1(t)
%           \  |
%            \ v
%  Qd(t-1) -> Qd(t)
%            /
%           /
%  Fd+1(t-1) 
%
% We assume parents are ordered (numbered) as follows:
% Qd(t-1), Fd+1(t-1), Fd(t-1), Q1(t), ..., Qd(t)
%
% The parents of Qd(t) can either be just Qd-1(t) or the whole stack Q1:d-1(t) (allQ)
% In either case, we will call them Qps.
% If d=1, Qps does not exist. Also, the F1(t-1) -> Q1(t) arc is optional.
% If the arc is missing, startprob does not need to be specified,
% since the toplevel is assumed to never reset (F1 does not exist).
% If d=D, Fd+1(t-1) does not exist (there is no signal from below).
%
% optional args [defaults]
%
% transprob - transprob(i,k,j) = prob transition from i to j given Qps = k ['leftright']
% selfprob  - prob of a transition from i to i given Qps=k [0.1]
% startprob - startprob(k,j) = prob start in j given Qps = k ['leftstart']
% startargs - other args to be passed to the sub tabular_CPD for learning startprob
% transargs - other args will be passed to the sub tabular_CPD for learning transprob
% allQ      - 1 means use all Q nodes above d as parents, 0 means just level d-1 [0]
% F1toQ1    - 1 means add F1(t-1) -> Q1(t) arc, 0 means level 1 never resets [0]
%
% For d=1, startprob(1,j) is only needed if F1toQ1=1
% Also, transprob(i,j) can be used instead of transprob(i,1,j).
%
% hhmmQ_CPD is a subclass of tabular_CPD so we inherit inference methods like CPD_to_pot, etc.
%
% We create isolated tabular_CPDs with no F parents to learn transprob/startprob
% so we can avail of e.g., entropic or Dirichlet priors.
% In the future, we will be able to represent the transprob using a tree_CPD.
%
% For details, see "Linear-time inference in hierarchical HMMs", Murphy and Paskin, NIPS'01.


ss = bnet.nnodes_per_slice;
%assert(self == Qnodes(d)+ss);
ns = bnet.node_sizes(:);
CPD.Qsizes = ns(Qnodes);
CPD.d = d;
CPD.D = D;
allQ = 0;

% find out which parents to use, to get right size
for i=1:2:length(varargin)
  switch varargin{i},
   case 'allQ', allQ = varargin{i+1}; 
  end
end

if d==1
  CPD.Qps = [];
else
  if allQ
    CPD.Qps = Qnodes(1:d-1);
  else
    CPD.Qps = Qnodes(d-1);
  end
end

Qsz = ns(self);
Qpsz = prod(ns(CPD.Qps));

% set default arguments
startprob = 'leftstart';
transprob = 'leftright';
startargs = {};
transargs = {};
CPD.F1toQ1 = 0;
selfprob = 0.1;

for i=1:2:length(varargin)
  switch varargin{i},
   case 'transprob', transprob = varargin{i+1}; 
   case 'selfprob',  selfprob = varargin{i+1}; 
   case 'startprob', startprob = varargin{i+1}; 
   case 'startargs', startargs = varargin{i+1}; 
   case 'transargs', transargs = varargin{i+1}; 
   case 'F1toQ1',    CPD.F1toQ1 = varargin{i+1}; 
  end
end

Qps = CPD.Qps + ss;
old_self = self-ss;

if strcmp(transprob, 'leftright')
  LR = mk_leftright_transmat(Qsz, selfprob);
  transprob = repmat(reshape(LR, [1 Qsz Qsz]), [Qpsz 1 1]); % transprob(k,i,j)
  transprob = permute(transprob, [2 1 3]); % now transprob(i,k,j)
end
transargs{end+1} = 'CPT';
transargs{end+1} = transprob;
CPD.sub_CPD_trans = mk_isolated_tabular_CPD([old_self Qps], ns([old_self Qps self]), transargs);
S = struct(CPD.sub_CPD_trans);
CPD.transprob = myreshape(S.CPT, [Qsz Qpsz Qsz]);


if strcmp(startprob, 'leftstart')
  startprob = zeros(Qpsz, Qsz);
  startprob(:,1) = 1;
end

if (d==1) & ~CPD.F1toQ1
  CPD.sub_CPD_start = [];
  CPD.startprob = [];
else
  startargs{end+1} = 'CPT';
  startargs{end+1} = startprob;
  CPD.sub_CPD_start = mk_isolated_tabular_CPD(Qps, ns([Qps self]), startargs);
  S = struct(CPD.sub_CPD_start);
  CPD.startprob = myreshape(S.CPT, [Qpsz Qsz]);
end

CPD = class(CPD, 'hhmmQ_CPD', tabular_CPD(bnet, self));

CPD = update_CPT(CPD);

