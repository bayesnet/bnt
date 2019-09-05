function engine = set_fields(engine, varargin)
% SET_FIELDS Set the fields for a generic engine
% engine = set_fields(engine, name/value pairs)
%
% e.g., engine = set_fields(engine, 'maximize', 1)

args = varargin;
nargs = length(args);
for i=1:2:nargs
  switch args{i},
   case 'maximize',
     engine.maximize = args{i+1};
     engine.jtree_engine = set_fields(engine.jtree_engine, 'maximize', args{i+1});
     engine.jtree_engine1 = set_fields(engine.jtree_engine1, 'maximize', args{i+1});
  end
end
