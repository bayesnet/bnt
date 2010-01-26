function [time, engine] = cmp_inference_dbn(bnet, engine, T, varargin)
% CMP_INFERENCE_DBN Compare several inference engines on a DBN
% function [time, engine] = cmp_inference_dbn(bnet, engine, T, ...)
%
% engine{i} is the i'th inference engine.
% time(e) = elapsed time for doing inference with engine e
%
% The list below gives optional arguments [default value in brackets].
%
% exact - specifies which engines do exact inference [ 1:length(engine) ]
% singletons_only - if 1, we only call marginal_nodes, else this  and marginal_family [0]
% check_ll - 1 means we check that the log-likelihoods are correct [1]

% set default params
exact = 1:length(engine);
singletons_only = 0;
check_ll = 1;
onodes = bnet.observed;

args = varargin;
nargs = length(args);
for i=1:2:nargs
  switch args{i},
   case 'exact', exact = args{i+1};
   case 'singletons_only', singletons_only = args{i+1};
   case 'check_ll', check_ll = args{i+1};
   case 'observed', onodes = args{i+1};
   otherwise,
    error(['unrecognized argument ' args{i}])
  end
end

E = length(engine);
ref = exact(1); % reference

ss = length(bnet.intra);
ev = sample_dbn(bnet, 'length', T);
evidence = cell(ss,T);
evidence(onodes,:) = ev(onodes, :);

for i=1:E
  tic;
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

if ~singletons_only
  get_marginals(engine, hnodes, exact, 0, T);
end
get_marginals(engine, hnodes, exact, 1, T);

%%%%%%%%%%

function get_marginals(engine, hnodes, exact, singletons, T)

bnet = bnet_from_engine(engine{1});
N = length(bnet.intra);
cnodes_bitv = zeros(1,N);
cnodes_bitv(bnet.cnodes) = 1;
ref = exact(1); % reference
cmp = exact(2:end);
E = length(engine);
m = cell(1,E);

for t=1:T
  for n=1:N
  %for n=hnodes(:)'
    for e=1:E
      if singletons
	m{e} = marginal_nodes(engine{e}, n, t);
      else
	m{e} = marginal_family(engine{e}, n, t);
      end
    end
    for e=cmp(:)'
      assert(isequal(m{e}.domain, m{ref}.domain));
      if cnodes_bitv(n) & isfield(m{e}, 'mu') & isfield(m{ref}, 'mu')
	wrong = ~approxeq(m{ref}.mu, m{e}.mu) | ~approxeq(m{ref}.Sigma, m{e}.Sigma);
      else
	wrong = ~approxeq(m{ref}.T(:), m{e}.T(:));
      end
      if wrong
	error(sprintf('engine %d is wrong; n=%d, t=%d, fam=%d', e, n, t, ~singletons))
      end
    end
  end
end
