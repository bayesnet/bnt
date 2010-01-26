function engine = jtree_unrolled_dbn_inf_engine(bnet, T, varargin)
% JTREE_UNROLLED_DBN_INF_ENGINE Unroll the DBN for T time-slices and apply jtree to the resulting static net
% engine = jtree_unrolled_dbn_inf_engine(bnet, T, ...)
%
% The following optional arguments can be specified in the form of name/value pairs:
% [default value in brackets]
%
% useC      - 1 means use jtree_C_inf_engine instead of jtree_inf_engine [0]
% constrained - 1 means we constrain ourselves to eliminate slice t before t+1 [1]
%
% e.g., engine = jtree_unrolled_inf_engine(bnet, 'useC', 1);

% set default params
N = length(bnet.intra);
useC = 0;
constrained = 1;

if nargin >= 3
  args = varargin;
  nargs = length(args);
  if isstr(args{1})
    for i=1:2:nargs
      switch args{i},
       case 'useC',   useC = args{i+1};
       case 'constrained',  constrained = args{i+1};
       otherwise,  
	error(['invalid argument name ' args{i}]);       
      end
    end
  else
    error(['invalid argument name ' args{1}]);       
  end
end

bnet2 = dbn_to_bnet(bnet, T);
ss = length(bnet.intra);
engine.ss = ss;

% If constrained_order = 1 we constrain ourselves to eliminate slice t before t+1.
% This prevents cliques containing nodes from far-apart time-slices.
if constrained
  stages = num2cell(unroll_set(1:ss, ss, T), 1);
else
  stages = { 1:length(bnet2.dag) };
end
if useC
  jengine = jtree_C_inf_engine(bnet2, 'stages', stages);
else
  jengine = jtree_inf_engine(bnet2, 'stages', stages);
end

engine.unrolled_engine = jengine;
% we don't inherit from jtree_inf_engine, because that would only store bnet2,
% and we would lose access to the DBN-specific fields like intra/inter

engine.nslices = T;
engine = class(engine, 'jtree_unrolled_dbn_inf_engine', inf_engine(bnet));
