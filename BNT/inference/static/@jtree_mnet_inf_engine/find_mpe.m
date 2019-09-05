function mpe = find_mpe(engine, evidence, varargin)
% FIND_MPE Find the most probable explanation of the data (assignment to the hidden nodes)
% function mpe = find_mpe(engine, evidence,...)
%
% evidence{i} = [] if X(i) is hidden, and otherwise contains its observed value (scalar or column vector).
%
% The following optional arguments can be specified in the form of name/value pairs:
% [default value in brackets]
%
% soft    - a cell array of soft/virtual evidence;
%           soft{i} is a prob. distrib. over i's values, or [] [ cell(1,N) ]
%

bnet = bnet_from_engine(engine);
ns = bnet.node_sizes(:);
N = length(bnet.dag);

engine.evidence = evidence;
  
% set default params
exclude = [];
soft_evidence = cell(1,N);

% parse optional params
args = varargin;
nargs = length(args);
for i=1:2:nargs
  switch args{i},
   case 'soft',    soft_evidence = args{i+1}; 
   otherwise,  
    error(['invalid argument name ' args{i}]);       
  end
end
engine.maximize = 1;

onodes = find(~isemptycell(evidence));
hnodes = find(isemptycell(evidence));
pot_type = determine_pot_type(bnet, onodes);
 if strcmp(pot_type, 'cg')
  check_for_cd_arcs(onodes, bnet.cnodes, bnet.dag);
end

hard_nodes = 1:N;
soft_nodes = find(~isemptycell(soft_evidence));
S = length(soft_nodes);
if S > 0
  assert(pot_type == 'd');
  assert(mysubset(soft_nodes, bnet.dnodes));
end
 
% Evaluate CPDs with evidence, and convert to potentials  
pot = cell(1, N+S);
for n=1:N
  fam = family(bnet.dag, n);
  e = bnet.equiv_class(n);
  if isempty(bnet.CPD{e})
    error(['must define CPD ' num2str(e)])
  else
    pot{n} = convert_to_pot(bnet.CPD{e}, pot_type, fam(:), evidence);
  end
end

for i=1:S
  n = soft_nodes(i);
  pot{N+i} = dpot(n, ns(n), soft_evidence{n});
end
clqs = engine.clq_ass_to_node([hard_nodes soft_nodes]); 

[clpot, seppot] = init_pot(engine, clqs, pot, pot_type, onodes);
[clpot, seppot] = collect_evidence(engine, clpot, seppot);
mpe = find_max_config(engine, clpot, seppot);
