function [engine, ll, niter] = enter_evidence(engine, evidence, varargin)
% ENTER_EVIDENCE Propagate evidence using belief propagation
% [engine, ll, niter] = enter_evidence(engine, evidence, ...)
%
% The log-likelihood is not computed; ll = 0.
% niter contains the number of iterations used (if engine.protocol = 'parallel')
%
% The following optional arguments can be specified in the form of name/value pairs:
% [default value in brackets]
%
% maximize - 1 means use max-product, 0 means use sum-product [0]
% exclude  - list of nodes whose potential will not be included in the joint [ [] ]
%
% e.g., engine = enter_evidence(engine, ev, 'maximize', 1)

ll = 0;
exclude = [];
maximize = 0;

if nargin >= 3
  args = varargin;
  nargs = length(args);
  for i=1:2:nargs
    switch args{i},
     case 'exclude', exclude = args{i+1};
     case 'maximize', maximize = args{i+1};
     otherwise,
      error(['invalid argument name ' args{i}]);
    end
  end
end

engine.maximize = maximize;

if ~isempty(engine.filename)
  engine.fid = fopen(engine.filename, 'w');
  if engine.fid == 0
    error(['can''t open ' engine.filename]);
  end
else
  engine.fid = [];
end

gdl = engine.gdl;
bnet = bnet_from_engine(engine);

ndoms = length(gdl.doms);
ns = bnet.node_sizes;
onodes = find(~isemptycell(evidence));
pot_type = determine_pot_type(bnet, onodes);

% prime each local kernel with evidence (if any)
local_kernel = cell(1, ndoms);
for i=1:ndoms
  if myismember(i, exclude)
    local_kernel{i} =  mk_initial_pot(pot_type, gdl.doms{i}, ns, bnet.cnodes, onodes);
  else
    e = bnet.equiv_class(i);
    local_kernel{i} =  convert_to_pot(bnet.CPD{e}, pot_type, gdl.doms{i}(:), evidence);
  end
end
  
% initialise all msgs to 1s
msg = cell(ndoms, ndoms);
for i=1:ndoms
  nbrs = gdl.nbrs{i};
  for j=nbrs(:)'
    dom = gdl.sepset{i,j};
    msg{i,j} = mk_initial_pot(pot_type, dom, ns, bnet.cnodes, onodes);
  end
end

switch engine.protocol
 case 'parallel', 
   [engine.marginal_domains, niter] = parallel_protocol(engine, evidence, pot_type, local_kernel, msg);
 case 'tree',
  engine.marginal_domains = serial_protocol(engine, evidence, pot_type, local_kernel, msg);
  niter = 1;
end
engine.niter = niter;

%fprintf('just finished %d iterations of belprop\n', niter);

if ~isempty(engine.filename)
  fclose(engine.fid);
end
