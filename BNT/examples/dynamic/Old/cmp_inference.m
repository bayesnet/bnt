function [err, time, engine] = cmp_inference(bnet, engine, exact, T, filter, singletons, maximize)
% CMP_INFERENCE Compare several inference engines on a DBN
% [err, time, engine] = cmp_inference(bnet, engine, exact, T, filter, singletons, maximize)
%
% engine{i} is the i'th inference engine.
% 'exact' specifies which engines do exact inference - 
%   we check that these all give the same results.
% 'T' is the length of the random sequence we generate.
% If filter=1, we do filtering, else smoothing (default: smoothing)
% If singletons=1, we compare marginal_nodes, else marginal_family (default: family)
%
% err(e,n,t) = sum_i | Pr_exact(X(n,t)=i) - Pr_e(X(n,t)=i) |
%   where Pr_e = prob. according to engine e
% time(e) = elapsed time for doing inference with engine e

err = [];

if nargin < 5, filter = 0; end
if nargin < 6, singletons = 0; end
if nargin < 7, maximize = 0; end

check_ll = 1;

assert(~maximize);

E = length(engine);
ref = exact(1); % reference

ss = length(bnet.intra);
ev = sample_dbn(bnet, 'length', T);
evidence = cell(ss,T);
onodes = bnet.observed;
evidence(onodes,:) = ev(onodes, :);

assert(~filter);
for i=1:E
  tic;
  %[engine{i}, ll(i)] = enter_evidence(engine{i}, evidence, 'maximize', maximize);
  [engine{i}, ll(i)] = enter_evidence(engine{i}, evidence);
  time(i)=toc;
  fprintf('engine %d took %6.4f seconds\n', i, time(i));
end

cmp = mysetdiff(exact, ref);
if check_ll
for i=cmp(:)'
  if ~approxeq(ll(ref), ll(i))
    error(['engine ' num2str(i) ' has wrong ll'])
  end
end
end
ll

hnodes = mysetdiff(1:ss, onodes);
m = cell(1,E);
for t=1:T
  for n=hnodes(:)'
    for e=1:E
      if singletons
	m{e} = marginal_nodes(engine{e}, n, t);
      else
	m{e} = marginal_family(engine{e}, n, t);
      end
    end
    for e=1:E
      assert(isequal(m{e}.domain, m{ref}.domain));
    end
    for e=cmp(:)'
      if ~approxeq(m{ref}.T(:), m{e}.T(:))
	str= sprintf('engine %d is wrong; n=%d, t=%d', e, n, t);
	error(str)
      end
    end
  end
end
