function [engine, loglik] = enter_evidence(engine, evidence, varargin)
% ENTER_EVIDENCE Add the specified evidence to the network (jtree_dbn)
% [engine, loglik] = enter_evidence(engine, evidence, ...)
%
% evidence{i,t} = [] if if X(i,t) is hidden, and otherwise contains its observed value (scalar or column vector)
%
% The following optional arguments can be specified in the form of name/value pairs:
% [default value in brackets]
%
% maximize - if 1, does max-product instead of sum-product [engine.maximize]
% softCPDpot{n,t} - use soft potential for node n instead of its CPD; set to [] to use CPD
% soft_evidence_nodes(i,1:2) = [n t] means the i'th piece of soft evidence is on node n in slice t 
% soft_evidence{i} - prob distribution over values for soft_evidence_nodes(i,:)
%
% e.g., engine = enter_evidence(engine, ev, 'maximize', 1)


% for add_ev in marginal_nodes
T = size(evidence, 2);
engine.evidence = evidence;
bnet = bnet_from_engine(engine);
ss = length(bnet.node_sizes_slice);
ns = bnet.node_sizes_slice(:);
engine.node_sizes = repmat(ns, [1 T]);
softCPDpot = cell(ss,T);
soft_evidence = {};
soft_evidence_nodes = [];

% parse optional params
args = varargin;
nargs = length(args);
if nargs > 0
  for i=1:2:nargs
    switch args{i},
     case 'maximize', engine.maximize = args{i+1}; 
     case 'softCPDpot', softCPDpot = args{i+1};
     case 'soft_evidence', soft_evidence = args{i+1};
     case 'soft_evidence_nodes', soft_evidence_nodes = args{i+1};
     otherwise,  
      error(['invalid argument name ' args{i}]);       
    end
  end
end

engine.jtree_engine = set_fields(engine.jtree_engine, 'maximize', engine.maximize);
engine.jtree_engine1 = set_fields(engine.jtree_engine1, 'maximize', engine.maximize);

[ss T] = size(evidence);
engine.T = T;
observed_bitv = ~isemptycell(evidence);
onodes = find(observed_bitv);
pot_type = determine_pot_type(bnet, onodes);
CPDpot = convert_dbn_CPDs_to_pots(bnet, evidence, pot_type, softCPDpot);

if ~isempty(soft_evidence_nodes)
  nsoft = size(soft_evidence_nodes,1);
  for i=1:nsoft
    n = soft_evidence_nodes(i,1);
    t = soft_evidence_nodes(i,2);
    if t==1
      dom = n;
    else
      dom = n+ss;
    end
    pot = dpot(dom, ns(n), soft_evidence{i});
    CPDpot{n,t} = multiply_by_pot(CPDpot{n,t}, pot);
  end
end
  
[engine.clpot, loglik] = enter_soft_evidence(engine, CPDpot, observed_bitv, pot_type);
