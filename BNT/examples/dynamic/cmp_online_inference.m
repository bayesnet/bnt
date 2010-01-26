function [time, engine] = cmp_online_inference(bnet, engine, T, varargin)
% CMP_ONLINE_INFERENCE Compare several online inference engines on a DBN
% function [time, engine] = cmp_online_inference(bnet, engine, T, ...)
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
cmp = mysetdiff(exact, ref);

ss = length(bnet.intra);
hnodes = mysetdiff(1:ss, onodes);
ev = sample_dbn(bnet, 'length', T);
evidence = cell(ss,T);
evidence(onodes,:) = ev(onodes, :);

time = zeros(1,E);
for t=1:T
  for e=1:E
    tic;
    [engine{e}, ll(e)] = enter_evidence(engine{e}, evidence(:,t), t);
    time(e)= time(e) + toc;
  end
  if check_ll
    for e=cmp(:)'
      if ~approxeq(ll(ref), ll(e))
	error(['engine ' num2str(e) ' has wrong ll'])
      end
    end
  end
  if ~singletons_only
    check_marginals(engine, hnodes, exact, 0, t);
  end
  check_marginals(engine, hnodes, exact, 1, t);
end


%%%%%%%%%%

function check_marginals(engine, hnodes, exact, singletons, t)

bnet = bnet_from_engine(engine{1});
N = length(bnet.intra);
cnodes_bitv = zeros(1,N);
cnodes_bitv(bnet.cnodes) = 1;
ref = exact(1); % reference
cmp = exact(2:end);
E = length(engine);
m = cell(1,E);

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
