function engine = likelihood_weighting_inf_engine(bnet, varargin)
% LIKELIHOOD_WEIGHTING_INF_ENGINE 
% engine = likelihood_weighting_inf_engine(bnet, ...)
%
% Optional arguments [defaults]
% nsamples - [500]

nsamples = 500;

if nargin >= 2
  args = varargin;
  nargs = length(args);
  for i=1:2:nargs
    switch args{i},
     case 'nsamples', nsamples= args{i+1};
     otherwise,
      error(['invalid argument name ' args{i}]);
    end
  end
end   

engine.nsamples = nsamples;
engine.samples = [];
engine.weights = [];
engine = class(engine, 'likelihood_weighting_inf_engine', inf_engine(bnet));
