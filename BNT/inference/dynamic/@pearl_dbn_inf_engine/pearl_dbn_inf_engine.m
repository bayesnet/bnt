function engine = pearl_dbn_inf_engine(bnet, varargin)
% LOOPY_DBN_INF_ENGINE Loopy Pearl version of forwards-backwards
% engine = loopy_dbn_inf_engine(bnet, ...)
%
% Optional arguments
% 'max_iter' - specifies the max num. forward-backward passes to perform [1]
% 'tol' - as in loopy_pearl [1e-3]
% 'momentum' - as in loopy_pearl [0]

error('pearl_dbn does not work yet')

max_iter = 1;
tol = 1e-3;
momentum = 0;

if nargin >= 2
  args = varargin;
  nargs = length(args);
  for i=1:2:nargs
    switch args{i},
     case 'max_iter', max_iter = args{i+1};
     case 'tol', tol = args{i+1};
     case 'momentum', momentum = args{i+1};
    end
  end
end

          
engine.max_iter = max_iter;
engine.tol = tol;
engine.momentum = momentum;
engine.pearl_engine = [];
engine.T = [];
engine.ss = length(bnet.intra);

engine.marginal = [];
engine.evidence = [];
engine.msg = [];
engine.parent_index = [];
engine.child_index = [];
%[engine.parent_index, engine.child_index] = mk_pearl_msg_indices(bnet); % need to unroll first

ss = length(bnet.intra);
engines.ss = ss;
onodes = bnet.observed;
hnodes = mysetdiff(1:ss, onodes);
obschild = zeros(1,ss);
for i=hnodes(:)'
  %ocs = myintersect(children(bnet.dag, i), onodes);
  ocs = children(bnet.intra, i);
  assert(length(ocs) <= 1);
  if length(ocs)==1
    obschild(i) = ocs(1);
  end
end
engine.obschild = obschild;

engine.mult_self_ndx = [];
engine.mult_parent_ndx = [];
engine.marg_self_ndx = [];
engine.marg_parent_ndx = [];


engine = class(engine, 'loopy_dbn_inf_engine', inf_engine(bnet));

