function [time, CPD, LL, cases] = cmp_learning_dbn(bnet, engine, T, varargin)
% CMP_LEARNING_DBN Compare a bunch of inference engines by learning a DBN
% function [time, CPD, LL, cases] = cmp_learning_dbn(bnet, engine, exact, T, ncases, max_iter)
%
% engine{i} is the i'th inference engine.
% time(e) = elapsed time for doing inference with engine e
% CPD{e,c} is the learned CPD for eclass c in engine e
% LL{e} is the learning curve for engine e
% cases{i} is the i'th training case
%
% The list below gives optional arguments [default value in brackets].
%
% exact - specifies which engines do exact inference [ 1:length(engine) ]
% check_ll - 1 means we check that the log-likelihoods are correct [1]
% ncases - num. random training cases [2]
% max_iter - max. num EM iterations [2]

% set default params
exact = 1:length(engine);
check_ll = 1;
ncases = 2;
max_iter = 2;

args = varargin;
nargs = length(args);
for i=1:2:nargs
  switch args{i},
   case 'exact', exact = args{i+1};
   case 'check_ll', check_ll = args{i+1};
   case 'ncases', ncases = args{i+1};
   case 'max_iter', max_iter = args{i+1};
   otherwise,
    error(['unrecognized argument ' args{i}])
  end
end

E = length(engine);
ss = length(bnet.intra);
onodes = bnet.observed;

cases = cell(1, ncases);
for i=1:ncases
  ev = sample_dbn(bnet, 'length', T);
  cases{i} = cell(ss,T);
  cases{i}(onodes,:) = ev(onodes, :);
end

LL = cell(1,E);
time = zeros(1,E);
for i=1:E
  tic
  [bnet2{i}, LL{i}] = learn_params_dbn_em(engine{i}, cases, 'max_iter', max_iter);
  time(i) = toc;
  fprintf('engine %d took %6.4f seconds\n', i, time(i));
end

ref = exact(1); % reference
cmp = mysetdiff(exact, ref);
if check_ll
  for i=cmp(:)'
    if ~approxeq(LL{ref}, LL{i})
      error(['engine ' num2str(i) ' has wrong ll'])
    end
  end
end

nCPDs = length(bnet.CPD);
CPD = cell(E, nCPDs);
tabular = zeros(1, nCPDs);
for i=1:E
  temp = bnet2{i};
  for c=1:nCPDs
    tabular(c) = isa(temp.CPD{c}, 'tabular_CPD');
    CPD{i,c} = struct(temp.CPD{c});
  end
end

for i=cmp(:)'
  for c=1:nCPDs
    if tabular(c)
      assert(approxeq(CPD{i,c}.CPT, CPD{ref,c}.CPT));
    else
      assert(approxeq(CPD{i,c}.mean, CPD{ref,c}.mean));
      assert(approxeq(CPD{i,c}.cov, CPD{ref,c}.cov));
      assert(approxeq(CPD{i,c}.weights, CPD{ref,c}.weights));
    end
  end
end

