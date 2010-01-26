function CPD = boolean_CPD(bnet, self, ftype, fname, pfail)
% BOOLEAN_CPD Make a tabular CPD representing a (noisy) boolean function
%
% CPD = boolean_cpd(bnet, self, 'inline', f) uses the inline function f
% to specify the CPT.
% e.g., suppose X4 = X2 AND (NOT X3). Then we can write
%    bnet.CPD{4} = boolean_CPD(bnet, 4, 'inline', inline('(x(1) & ~x(2)'));  
% Note that x(1) refers pvals(1) = X2, and x(2) refers to pvals(2)=X3.
%
% CPD = boolean_cpd(bnet, self, 'named', f) assumes f is a function name.
% f can be built-in to matlab, or a file.
% e.g., If X4 = X2 AND X3, we can write
%    bnet.CPD{4} = boolean_CPD(bnet, 4, 'named', 'and');
% e.g., If X4 = X2 OR X3, we can write
%    bnet.CPD{4} = boolean_CPD(bnet, 4, 'named', 'any');
%
% CPD = boolean_cpd(bnet, self, 'rnd') makes a random non-redundant bool fn.
%
% CPD = boolean_CPD(bnet, self, 'inline'/'named', f, pfail)
% will put probability mass 1-pfail on f(parents), and put pfail on the other value.
% This is useful for simulating noisy boolean functions.
% If pfail is omitted, it is set to 0.
% (Note that adding noise to a random (non-redundant) boolean function just creates a different
% (potentially redundant) random boolean function.)
%
% Note: This cannot be used to simulate a noisy-OR gate.
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
% By contrast, boolean_CPD(bnet, C, 'any', p) would define
%
%  A  B  P(C=0) 
%  0  0  1-p    
%  1  0  p      
%  0  1  p
%  1  1  p


if nargin==0
  % This occurs if we are trying to load an object from a file.
  CPD = tabular_CPD(bnet, self);
  return;
elseif isa(bnet, 'boolean_CPD')
  % This might occur if we are copying an object.
  CPD = bnet;
  return;
end

if nargin < 5, pfail = 0; end

ps = parents(bnet.dag, self);
ns = bnet.node_sizes;
psizes = ns(ps);
self_size = ns(self);

psucc = 1-pfail;

k = length(ps);
switch ftype
 case 'inline', f = eval_bool_fn(fname, k);
 case 'named',  f = eval_bool_fn(fname, k);
 case 'rnd',    f = mk_rnd_bool_fn(k);
 otherwise,     error(['unknown function type ' ftype]);
end

CPT = zeros(prod(psizes), self_size);
ndx = find(f==0);
CPT(ndx, 1) = psucc;
CPT(ndx, 2) = pfail;
ndx = find(f==1);
CPT(ndx, 2) = psucc;
CPT(ndx, 1) = pfail;
if k > 0
  CPT = reshape(CPT, [psizes self_size]);  
end

clamp = 1;
CPD = tabular_CPD(bnet, self, CPT, [], clamp);



%%%%%%%%%%%%

function f = eval_bool_fn(fname, n)
% EVAL_BOOL_FN Evaluate a boolean function on all bit vectors of length n
% f = eval_bool_fn(fname, n)
%
% e.g. f = eval_bool_fn(inline('x(1) & x(3)'), 3)
% returns   0     0     0     0     0     1     0     1

ns = 2*ones(1, n);
f = zeros(1, 2^n);
bits = ind2subv(ns, 1:2^n);
for i=1:2^n
  f(i) = feval(fname, bits(i,:)-1);
end

%%%%%%%%%%%%%%%

function f = mk_rnd_bool_fn(n)
% MK_RND_BOOL_FN Make a random bit vector of length n that encodes a non-redundant boolean function
% f = mk_rnd_bool_fn(n)

red = 1;
while red
  f = sample_discrete([0.5 0.5], 2^n, 1)-1;
  red = redundant_bool_fn(f);
end

%%%%%%%%


function red = redundant_bool_fn(f)
% REDUNDANT_BOOL_FN Does a boolean function depend on all its input values?
% r = redundant_bool_fn(f)
%
% f is a vector of length 2^n, representing the output for each bit vector.
% An input is redundant if there is no assignment to the other bits
% which changes the output e.g., input 1 is redundant if u(2:n) s.t.,
% f([0 u(2:n)]) <> f([1 u(2:n)]). 
% A function is redundant it it has any redundant inputs.

n = log2(length(f));
ns = 2*ones(1,n);
red = 0;
for i=1:n
  ens = ns;
  ens(i) = 1;
  U = ind2subv(ens, 1:2^(n-1));
  U(:,i) = 1;
  f1 = f(subv2ind(ns, U));
  U(:,i) = 2;
  f2 = f(subv2ind(ns, U));
  if isequal(f1, f2)
    red = 1;
    return;
  end
end


%%%%%%%%%%

function [b, iter] = rnd_truth_table(N)
% RND_TRUTH_TABLE Construct the output of a random truth table s.t. each input is non-redundant
% b = rnd_truth_table(N)
%
% N is the number of inputs. 
% b is a random bit string of length N, representing the output of the truth table.
% Non-redundant means that, for each input position k,
% there are at least two bit patterns, u and v, that differ only in the k'th position,
% s.t., f(u) ~= f(v), where f is the function represented by b.
% We use rejection sampling to ensure non-redundancy.
%
% Example: b = [0 0 0 1  0 0 0 1] is indep of 3rd input (AND of inputs 1 and 2)

bits = ind2subv(2*ones(1,N), 1:2^N)-1;
redundant = 1;
iter = 0;
while redundant & (iter < 4)
  iter = iter + 1;
  b = sample_discrete([0.5 0.5], 1, 2^N)-1;
  redundant = 0;
  for i=1:N
    on = find(bits(:,i)==1);
    off = find(bits(:,i)==0);
    if isequal(b(on), b(off))
      redundant = 1;
      break;
    end
  end
end

