function [time, engine] = cmp_inference_static(bnet, engine, varargin)
% CMP_INFERENCE Compare several inference engines on a BN
% function [time, engine] = cmp_inference_static(bnet, engine, ...)
%
% engine{i} is the i'th inference engine.
% time(e) = elapsed time for doing inference with engine e
%
% The list below gives optional arguments [default value in brackets].
%
% exact - specifies which engines do exact inference [ 1:length(engine) ]
% singletons_only - if 1, we only call marginal_nodes, else this  and marginal_family [0]
% maximize - 1 means we do max-propagation, 0 means sum-propagation [0]
% check_ll - 1 means we check that the log-likelihoods are correct [1]
% observed - list of the observed ndoes [ bnet.observed ]
% check_converged - list of loopy engines that should be checked for convergence [ [] ]
%    If an engine has converged, it is added to the exact list.


% set default params
exact = 1:length(engine);
singletons_only = 0;
maximize = 0;
check_ll = 1;
observed = bnet.observed;
check_converged = [];

args = varargin;
nargs = length(args);
for i=1:2:nargs
  switch args{i},
   case 'exact', exact = args{i+1};
   case 'singletons_only', singletons_only = args{i+1};
   case 'maximize', maximize = args{i+1};
   case 'check_ll', check_ll = args{i+1};
   case 'observed', observed = args{i+1};
   case 'check_converged', check_converged = args{i+1};
   otherwise,
    error(['unrecognized argument ' args{i}])
  end
end

E = length(engine);
ref = exact(1); % reference

N = length(bnet.dag);
ev = sample_bnet(bnet);
evidence = cell(1,N);
evidence(observed) = ev(observed);
%celldisp(evidence(observed))

for i=1:E
  tic;
  if check_ll
    [engine{i}, ll(i)] = enter_evidence(engine{i}, evidence, 'maximize', maximize);
  else
    engine{i} = enter_evidence(engine{i}, evidence, 'maximize', maximize);
  end
  time(i)=toc;
end

for i=check_converged(:)'
  niter = loopy_converged(engine{i});
  if niter > 0
    fprintf('loopy engine %d  converged in %d iterations\n', i, niter);
%    exact = myunion(exact, i);
  else
    fprintf('loopy engine %d has not converged\n', i);
  end
end

cmp = exact(2:end);
if check_ll
  for i=cmp(:)'
    assert(approxeq(ll(ref), ll(i)));
  end
end

hnodes = mysetdiff(1:N, observed);

if ~singletons_only
  get_marginals(engine, hnodes, exact, 0);
end
get_marginals(engine, hnodes, exact, 1);

%%%%%%%%%%

function get_marginals(engine, hnodes, exact, singletons)

bnet = bnet_from_engine(engine{1});
N = length(bnet.dag);
cnodes_bitv = zeros(1,N);
cnodes_bitv(bnet.cnodes) = 1;
ref = exact(1); % reference
cmp = exact(2:end);
E = length(engine);

for n=hnodes(:)'
  for e=1:E
    if singletons
      m{e} = marginal_nodes(engine{e}, n);
    else
      m{e} = marginal_family(engine{e}, n);
    end
  end
  for e=cmp(:)'
    if cnodes_bitv(n)
      assert(approxeq(m{ref}.mu, m{e}.mu))
      assert(approxeq(m{ref}.Sigma, m{e}.Sigma))
    else
      assert(approxeq(m{ref}.T, m{e}.T))
    end
    assert(isequal(m{e}.domain, m{ref}.domain));
  end
end
