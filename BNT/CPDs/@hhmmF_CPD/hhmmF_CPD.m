function CPD = hhmmF_CPD(bnet, self, Qself, Fbelow, varargin)
% HHMMF_CPD Make the CPD for an F node in a hierarchical HMM
% CPD = hhmmF_CPD(bnet, self, Qself,  Fbelow, ...)
%
%        Qps
%          \
%           \
%           Fself
%         /   |
%        /    |
%       Qself Fbelow
%
% We assume nodes are ordered (numbered) as follows: Qps, Q, Fbelow, F
% All nodes numbers should be from slice 1.
%
% If Fbelow if missing, this becomes a regular tabular_CPD.
% Qps may be omitted.
%
% optional args [defaults]
% 
% Qps - node numbers.
% termprob - termprob(k,i,2) = prob finishing given Q(d)=i and Q(1:d-1)=k [ finish in last state wp 0.9]
%
% hhmmF_CPD is a subclass of tabular_CPD so we inherit inference methods like CPD_to_pot, etc.
%
% We create an isolated tabular_CPD with no F parent to learn termprob
% so we can avail of e.g., entropic or Dirichlet priors.
%
% For details, see "Linear-time inference in hierarchical HMMs", Murphy and Paskin, NIPS'01.



Qps = [];
% get parents
for i=1:2:length(varargin)
  switch varargin{i},
   case 'Qps', Qps = varargin{i+1}; 
  end
end

ns = bnet.node_sizes(:);
Qsz = ns(Qself);
Qpsz = prod(ns(Qps));
CPD.Qsz = Qsz;
CPD.Qpsz = Qpsz;

ps = parents(bnet.dag, self);
CPD.Fbelow_ndx = find_equiv_posns(Fbelow, ps);
CPD.Qps_ndx = find_equiv_posns(Qps, ps);
CPD.Qself_ndx = find_equiv_posns(Qself, ps);

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
  end
end

CPD.sub_CPD_term = mk_isolated_tabular_CPD([Qpsz Qsz 2], {'CPT', termprob});
S = struct(CPD.sub_CPD_term);
CPD.termprob = S.CPT;

CPD = class(CPD, 'hhmmF_CPD', tabular_CPD(bnet, self));

CPD = update_CPT(CPD);

