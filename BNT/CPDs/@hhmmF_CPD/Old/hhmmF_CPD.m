function CPD = hhmmF_CPD(bnet, self, Qnodes, d, D, varargin)
% HHMMF_CPD Make the CPD for an F node at depth D of a D-level hierarchical HMM
% CPD = hhmmF_CPD(bnet, self, Qnodes, d, D, ...)
%
%    Q(d-1)
%          \
%           \
%           F(d)
%         /   |
%        /    |
%    Q(d)  F(d+1)
%
% We assume nodes are ordered (numbered) as follows:
% Q(1), ... Q(d), F(d+1), F(d)
%
% F(d)=2 means level d has finished. The prob this happens depends on Q(d)
% and optionally on Q(d-1), Q(d=1), ..., Q(1).
% Also, level d can only finish if the level below has finished
% (hence the F(d+1) -> F(d) arc).
%
% If d=D, there is no F(d+1), so F(d) is just a regular tabular_CPD.
% If all models always finish in the same state (e.g., their last),
% we don't need to condition on the state of parent models (Q(d-1), ...)
%
% optional args [defaults]
%
% termprob - termprob(k,i,2) = prob finishing given Q(d)=i and Q(1:d-1)=k [ finish in last state ]
%
% hhmmF_CPD is a subclass of tabular_CPD so we inherit inference methods like CPD_to_pot, etc.
%
% We create an isolated tabular_CPD with no F parent to learn termprob
% so we can avail of e.g., entropic or Dirichlet priors.
%
% For details, see "Linear-time inference in hierarchical HMMs", Murphy and Paskin, NIPS'01.


ps = parents(bnet.dag, self);
Qps = myintersect(ps, Qnodes);
F = mysetdiff(ps, Qps);
CPD.Q = Qps(end); % Q(d)
assert(CPD.Q == Qnodes(d));
CPD.Qps = Qps(1:end-1); % all Q parents except Q(d), i.e., calling context

ns = bnet.node_sizes(:);
CPD.Qsizes = ns(Qnodes);
CPD.d = d;
CPD.D = D;

Qsz = ns(CPD.Q);
Qpsz = prod(ns(CPD.Qps));

% set default arguments
p = 0.9;
%termprob(k,i,t) Might terminate if i=Qsz; will not terminate if i<Qsz
termprob = zeros(Qpsz, Qsz, 2);
termprob(:, Qsz, 2) = p; 
termprob(:, Qsz, 1) = 1-p; 
termprob(:, 1:(Qsz-1), 1) = 1; 
    
for i=1:2:length(varargin)
  switch varargin{i},
   case 'termprob', termprob = varargin{i+1}; 
   otherwise, error(['unrecognized argument ' varargin{i}])
  end
end

ps = [CPD.Qps CPD.Q];
% ns(self) = 2 since this is an F node
CPD.sub_CPD_term = mk_isolated_tabular_CPD(ps, ns([ps self]), {'CPT', termprob});
S = struct(CPD.sub_CPD_term);
CPD.termprob = S.CPT;

CPD = class(CPD, 'hhmmF_CPD', tabular_CPD(bnet, self));

CPD = update_CPT(CPD);

