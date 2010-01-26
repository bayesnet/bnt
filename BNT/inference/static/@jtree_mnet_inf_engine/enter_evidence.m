function [engine, loglik] = enter_evidence(engine, evidence, varargin)
% ENTER_EVIDENCE Add the specified evidence to the network (jtree)
% [engine, loglik] = enter_evidence(engine, evidence, ...)
%
% evidence{i} = [] if X(i) is hidden, and otherwise contains its observed value (scalar or column vector).
%
% The following optional arguments can be specified in the form of name/value pairs:
% [default value in brackets]
%
% soft    - a cell array of soft/virtual evidence;
%           soft{i} is a prob. distrib. over i's values, or [] [ cell(1,N) ]
%
% e.g., engine = enter_evidence(engine, ev, 'soft', soft_ev)

bnet = bnet_from_engine(engine);
ns = bnet.node_sizes(:);
N = length(bnet.dag);

engine.evidence = evidence; % store this for marginal_nodes with add_ev option
engine.maximize = 0;

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

onodes = find(~isemptycell(evidence));
hnodes = find(isemptycell(evidence));
pot_type = determine_pot_type(bnet, onodes);
 if strcmp(pot_type, 'cg')
  check_for_cd_arcs(onodes, bnet.cnodes, bnet.dag);
end

% Evaluate CPDs with evidence, and convert to potentials  
pot = cell(1, N);
for n=1:N
  fam = family(bnet.dag, n);
  e = bnet.equiv_class(n);
  if isempty(bnet.CPD{e})
    error(['must define CPD ' num2str(e)])
  else
    pot{n} = convert_to_pot(bnet.CPD{e}, pot_type, fam(:), evidence);
  end
end
clqs = engine.clq_ass_to_node(1:N);

% soft evidence
soft_nodes = find(~isemptycell(soft_evidence));
S = length(soft_nodes);
if S > 0
  assert(pot_type == 'd');
  assert(mysubset(soft_nodes, bnet.dnodes));
end
for i=1:S
  n = soft_nodes(i);
  pot{end+1} = dpot(n, ns(n), soft_evidence{n});
end
clqs = [clqs engine.clq_ass_to_node(soft_nodes)]; 


[clpot, seppot] = init_pot(engine, clqs, pot, pot_type, onodes);
[clpot, seppot] = collect_evidence(engine, clpot, seppot);
[clpot, seppot] = distribute_evidence(engine, clpot, seppot);

C = length(clpot);
ll = zeros(1, C);
for i=1:C
  [clpot{i}, ll(i)] = normalize_pot(clpot{i});
end
loglik = ll(1); % we can extract the likelihood from any clique

engine.clpot = clpot;
