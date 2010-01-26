function CPD = noisyor_CPD(bnet, self, leak_inhibit, inhibit)
% NOISYOR_CPD Make a noisy-or CPD
% CPD = NOISYOR_CPD(BNET, NODE_NUM, LEAK_INHIBIT, INHIBIT)
%
% A noisy-or node turns on if any of its parents are on, provided they are not inhibited.
% The prob. that the i'th parent gets inhibited (flipped from 1 to 0) is inhibit(i).
% The prob that the leak node (a dummy parent that is always on) gets inhibit is leak_inhibit.
% These params default to random values if omitted.
%
% Example: suppose C has parents A and B, and the
% link of A->C fails with prob pA and the link B->C fails with pB.
% Then the noisy-OR gate defines the following distribution
%
%  A  B  P(C=0)
%  0  0  1.0
%  1  0  pA
%  0  1  pB
%  1  1  pA * PB
%
% Currently, learning is not supported for noisy-or nodes
% (since the M step is somewhat complicated).
%
% For simple generalizations of the noisy-OR model, see e.g.,
% - Srinivas, "A generalization of the noisy-OR model", UAI 93
% - Meek and Heckerman, "Learning Causal interaction models", UAI 97.
  


if nargin==0
  % This occurs if we are trying to load an object from a file.
  CPD = init_fields;
  CPD = class(CPD, 'noisyor_CPD', discrete_CPD(1, []));
  return;
elseif isa(bnet, 'noisyor_CPD')
  % This might occur if we are copying an object.
  CPD = bnet;
  return;
end
CPD = init_fields;


ps = parents(bnet.dag, self);
fam = [ps self];
ns = bnet.node_sizes;
assert(all(ns(fam)==2));
assert(isempty(myintersect(fam, bnet.cnodes)));

if nargin < 3, leak_inhibit = rand(1, 1); end
if nargin < 4, inhibit = rand(1, length(ps)); end

CPD.self = self;
CPD.inhibit = inhibit;
CPD.leak_inhibit = leak_inhibit;


% For BIC
CPD.nparams = 0;
CPD.nsamples = 0;

CPD.CPT = []; % cached copy, to speed up CPD_to_CPT

clamped = 1;
CPD = class(CPD, 'noisyor_CPD', discrete_CPD(clamped, ns([ps self])));



%%%%%%%%%%%

function CPD = init_fields()
% This ensures we define the fields in the same order 
% no matter whether we load an object from a file,
% or create it from scratch. (Matlab requires this.)

CPD.self = [];
CPD.inhibit = [];
CPD.leak_inhibit = [];
CPD.nparams = [];
CPD.nsamples = [];
CPD.CPT = [];
