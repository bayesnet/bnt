function CPD = hhmmQ_CPD(bnet, self, varargin)
% HHMMQ_CPD Make the CPD for a Q node in a hierarchical HMM
% CPD = hhmmQ_CPD(bnet, self, ...)
%
%  Fself(t-1)   Qps(t)
%           \    |
%            \   v
%  Qold(t-1) ->  Q(t)
%            /
%           /
%  Fbelow(t-1) 
%
% Let ss = slice size = num. nodes per slice.
% This node is Q(t), and has mandatory parents Qold(t-1) (assumed to be numbered Q(t)-ss)
% and optional parents Fbelow, Fself, Qps.
% We require parents to be ordered (numbered) as follows:
% Qold, Fbelow, Fself, Qps, Q.
%
% If Fself=2, we use the transition matrix, else we use the prior matrix.
% If Fself node is omitted (eg. top level), we always use the transition matrix.
% If Fbelow=2, we may change state, otherwise we must stay in the same state.
% If Fbelow node is omitted (eg., bottom level), we may change state at every step.
% If Qps (Q parents) are specified, all parameters are conditioned on their joint value.
% We may choose any subset of nodes to condition on, as long as they as numbered lower than self.
%
% optional args [defaults]
%
% Fself - node number <= ss
% Fbelow  - node number  <= ss
% Qps - node numbers (all <= 2*ss) - uses 2TBN indexing
% transprob - transprob(i,k,j) = prob transition from i to j given Qps = k ['leftright']
% selfprob  - prob of a transition from i to i given Qps=k [0.1]
% startprob - startprob(k,j) = prob start in j given Qps = k ['leftstart']
% startargs - other args to be passed to the sub tabular_CPD for learning startprob
% transargs - other args will be passed to the sub tabular_CPD for learning transprob
% fullstartprob - 1 means startprob depends on Q(t-1) [0]
% hhmmQ_CPD is a subclass of tabular_CPD so we inherit inference methods like CPD_to_pot, etc.
%
% We create isolated tabular_CPDs with no F parents to learn transprob/startprob
% so we can avail of e.g., entropic or Dirichlet priors.
% In the future, we will be able to represent the transprob using a tree_CPD.
%
% For details, see "Linear-time inference in hierarchical HMMs", Murphy and Paskin, NIPS'01.


ss = bnet.nnodes_per_slice;
ns = bnet.node_sizes(:);

% set default arguments
Fself = [];
Fbelow = [];
Qps = [];
startprob = 'leftstart';
transprob = 'leftright';
startargs = {};
transargs = {};
selfprob = 0.1;
fullstartprob = 0;

for i=1:2:length(varargin)
  switch varargin{i},
   case 'Fself', Fself = varargin{i+1};
   case 'Fbelow', Fbelow = varargin{i+1};
   case 'Qps', Qps = varargin{i+1};
   case 'transprob', transprob = varargin{i+1}; 
   case 'selfprob',  selfprob = varargin{i+1}; 
   case 'startprob', startprob = varargin{i+1}; 
   case 'startargs', startargs = varargin{i+1}; 
   case 'transargs', transargs = varargin{i+1}; 
   case 'fullstartprob', fullstartprob = varargin{i+1}; 
  end
end

CPD.fullstartprob = fullstartprob;

ps = parents(bnet.dag, self);
ndsz = ns(:)';
CPD.dom_sz = [ndsz(ps) ns(self)];
CPD.Fself_ndx = find_equiv_posns(Fself, ps);
CPD.Fbelow_ndx = find_equiv_posns(Fbelow, ps);
%CPD.Qps_ndx = find_equiv_posns(Qps+ss, ps);
CPD.Qps_ndx = find_equiv_posns(Qps, ps);
old_self = self-ss;
CPD.old_self_ndx = find_equiv_posns(old_self, ps);

Qps = ps(CPD.Qps_ndx);
CPD.Qsz = ns(self);
CPD.Qpsz = prod(ns(Qps));
CPD.Qpsizes = ns(Qps);
Qsz = CPD.Qsz;
Qpsz = CPD.Qpsz;

if strcmp(transprob, 'leftright')
  LR = mk_leftright_transmat(Qsz, selfprob);
  transprob = repmat(reshape(LR, [1 Qsz Qsz]), [Qpsz 1 1]); % transprob(k,i,j)
  transprob = permute(transprob, [2 1 3]); % now transprob(i,k,j)
end
transargs{end+1} = 'CPT';
transargs{end+1} = transprob;
CPD.sub_CPD_trans = mk_isolated_tabular_CPD(ns([old_self Qps self]), transargs);
S = struct(CPD.sub_CPD_trans);
%CPD.transprob = myreshape(S.CPT, [Qsz Qpsz Qsz]);
CPD.transprob = S.CPT;


if strcmp(startprob, 'leftstart')
  startprob = zeros(Qpsz, Qsz);
  startprob(:,1) = 1;
end
if isempty(CPD.Fself_ndx)
  CPD.sub_CPD_start = [];
  CPD.startprob = [];
else
  startargs{end+1} = 'CPT';
  startargs{end+1} = startprob;
  if CPD.fullstartprob
    CPD.sub_CPD_start = mk_isolated_tabular_CPD(ns([self Qps self]), startargs);
    S = struct(CPD.sub_CPD_start);
    %CPD.startprob = myreshape(S.CPT, [Qsz Qpsz Qsz]);
    CPD.startprob = S.CPT;
  else
    CPD.sub_CPD_start = mk_isolated_tabular_CPD(ns([Qps self]), startargs);
    S = struct(CPD.sub_CPD_start);
    %CPD.startprob = myreshape(S.CPT, [CPD.Qpsizes Qsz]);
    CPD.startprob = S.CPT;
  end
end

CPD = class(CPD, 'hhmmQ_CPD', tabular_CPD(bnet, self));

CPD = update_CPT(CPD);

