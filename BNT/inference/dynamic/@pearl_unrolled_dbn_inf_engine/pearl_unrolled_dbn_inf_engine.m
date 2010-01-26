function engine = pearl_unrolled_dbn_inf_engine(bnet, varargin)
% LOOPY_DBN_INF_ENGINE Loopy Pearl version of forwards-backwards
% engine = loopy_unrolld_dbn_inf_engine(bnet, ...)
%
% Optional arguments
% 'max_iter' - specifies the max num. forward-backward passes to perform PER SLICE [2]
% 'tol' - as in loopy_pearl [1e-3]
% 'momentum' - as in loopy_pearl [0]
% protocol - tree or parallel [parallel]
% filename - as in pearl [ '' ]

max_iter_per_slice = 2;
tol = 1e-3;
momentum = 0;
protocol = 'parallel';
filename = '';

args = varargin;
for i=1:2:length(args)
  switch args{i},
   case 'max_iter', max_iter_per_slice = args{i+1};
   case 'tol', tol = args{i+1};
   case 'momentum', momentum = args{i+1};
   case 'protocol', protocol = args{i+1};
   case 'filename', filename = args{i+1};
  end
end

engine.filename = filename;
engine.max_iter_per_slice = max_iter_per_slice;
engine.tol = tol;
engine.momentum = momentum;
engine.unrolled_engine = [];
engine.T = -1;
engine.ss = length(bnet.intra);
engine.protocol = protocol;

engine = class(engine, 'pearl_unrolled_dbn_inf_engine', inf_engine(bnet));

