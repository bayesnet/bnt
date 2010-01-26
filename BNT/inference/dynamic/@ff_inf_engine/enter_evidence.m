function [engine, loglik] = enter_evidence(engine, evidence, varargin)
% ENTER_EVIDENCE Add the specified evidence to the network (ff)
% [engine, loglik] = enter_evidence(engine, evidence, ...)
%
% evidence{i,t} = [] if if X(i,t) is hidden, and otherwise contains its observed value (scalar or
% column vector)
%
% The following optional arguments can be specified in the form of name/value pairs:
% [default value in brackets]
%
% maximize - if 1, does max-product (not yet supported), else sum-product [0]
% filter -   if 1, do filtering, else smoothing [0]
%
% e.g., engine = enter_evidence(engine, ev, 'maximize', 1)

maximize = 0;
filter = 0;

% parse optional params
args = varargin;
nargs = length(args);
if nargs > 0
  for i=1:2:nargs
    switch args{i},
     case 'maximize', maximize = args{i+1}; 
     case 'filter', filter = args{i+1}; 
     otherwise,  
      error(['invalid argument name ' args{i}]);       
    end
  end
end

assert(~maximize);


[ss T] = size(evidence);
observed = ~isemptycell(evidence);
bnet = bnet_from_engine(engine);
%pot_type = determine_pot_type(find(observed(:,1)), bnet.cnodes_slice, bnet.intra);
pot_type = determine_pot_type(bnet, observed);
% we assume we can use the same pot_type in all slices

CPDpot = convert_dbn_CPDs_to_pots(bnet, evidence, pot_type);

% Now convert CPDs on observed nodes to be potentials just on their parents
assert(pot_type == 'd');
onodes = bnet.observed(:);
ns = bnet.node_sizes_slice;
ns(onodes) = 1;
for t=1:T
  for i=onodes
    p = parents(bnet.dag, i);
    %CPDpot{i,t} = set_domain_pot(CPDpot{i,t}, p); % leaves size too long
    temp = pot_to_marginal(CPDpot{i,t});
    CPDpot{i,t} = dpot(p, ns(p), temp.T); % assumes pot_type = d
  end
end

[engine.marginals, engine.fwd, engine.back, loglik] = enter_soft_evidence(engine, CPDpot, observed, pot_type, filter);

engine.CPDpot = CPDpot;
engine.filter = filter;
