function [engine, loglik] = enter_evidence(engine, evidence, varargin)
% ENTER_EVIDENCE Add the specified evidence to the network (jtree)
% [engine, loglik] = enter_evidence(engine, evidence, ...)
%
% evidence{i} = [] if X(i) is hidden, and otherwise contains its observed value (scalar or column vector).
%
% The following optional arguments can be specified in the form of name/value pairs:
% [default value in brackets]
%
% maximize - if 1, does max-product instead of sum-product [0]
% soft    - a cell array of soft/virtual evidence;
%           soft{i} is a prob. distrib. over i's values, or [] [ cell(1,N) ]
%
% e.g., engine = enter_evidence(engine, ev, 'soft', soft_ev)
%
% For backwards compatibility with BNT2, you can also specify the parameters in the following order
%  engine = enter_evidence(engine, ev, soft_ev)

bnet = bnet_from_engine(engine);
ns = bnet.node_sizes(:);
N = length(bnet.dag);

engine.evidence = evidence; % store this for marginal_nodes with add_ev option
  
% set default params
exclude = [];
soft_evidence = cell(1,N);
maximize = 0;

% parse optional params
args = varargin;
nargs = length(args);
if nargs > 0
  if iscell(args{1})
    soft_evidence = args{1};
  else
    for i=1:2:nargs
      switch args{i},
       case 'soft',    soft_evidence = args{i+1}; 
       case 'maximize', maximize = args{i+1}; 
       otherwise,  
	error(['invalid argument name ' args{i}]);       
      end
    end
  end
end

engine.maximize = maximize;

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
[clpot, seppot] = distribute_evidence(engine, clpot, seppot);
C = length(clpot);
ll = zeros(1, C);
for i=1:C
   domain = clpot{i}.domain;
   sizes = clpot{i}.sizes;
   T = clpot{i}.T;
   clpot{i} = dpot(domain, sizes, T);
end
   
for i=1:C
  [clpot{i}, ll(i)] = normalize_pot(clpot{i});
end
loglik = ll(1); % we can extract the likelihood from any clique

engine.clpot = clpot;
