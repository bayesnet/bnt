function [engine, ll, niter] = enter_evidence(engine, evidence, varargin)
% ENTER_EVIDENCE Propagate evidence using belief propagation
% [engine, ll, niter] = enter_evidence(engine, evidence, ...)
%
% The log-likelihood is not computed; ll = 0.
% niter contains the number of iterations used 
%
% The following optional arguments can be specified in the form of name/value pairs:
% [default value in brackets]
%
% maximize - 1 means use max-product, 0 means use sum-product [0]
%
% e.g., engine = enter_evidence(engine, ev, 'maximize', 1)

ll = 0;
maximize = 0;

if nargin >= 3
  args = varargin;
  nargs = length(args);
  for i=1:2:nargs
    switch args{i},
     case 'maximize', maximize = args{i+1};
     otherwise,
      error(['invalid argument name ' args{i}]);
    end
  end
end

verbose = 0;

ns = engine.fgraph.node_sizes;
onodes = find(~isemptycell(evidence));
hnodes = find(isemptycell(evidence));
cnodes = engine.fgraph.cnodes;
pot_type = determine_pot_type(engine.fgraph, onodes);

% prime each local kernel with evidence (if any)
nfactors = engine.fgraph.nfactors;
nvars = engine.fgraph.nvars;
factors = cell(1,nfactors);
for f=1:nfactors
  K = engine.fgraph.factors{engine.fgraph.equiv_class(f)};
  factors{f} = convert_to_pot(K, pot_type, engine.fgraph.dom{f}(:), evidence);
end
  
% initialise msgs
msg_var_to_fac = cell(nvars, nfactors);
for x=1:nvars
  for f=engine.fgraph.dep{x}
    msg_var_to_fac{x,f} = mk_initial_pot(pot_type, x, ns, cnodes, onodes);
  end
end
msg_fac_to_var = cell(nfactors, nvars);
dom = cell(1, nfactors);
for f=1:nfactors
  %hdom{f} = myintersect(engine.fgraph.dom{f}, hnodes);
  dom{f} = engine.fgraph.dom{f}(:)';
  for x=dom{f}
    msg_fac_to_var{f,x} = mk_initial_pot(pot_type, x, ns, cnodes, onodes);
    %msg_fac_to_var{f,x} = marginalize_pot(factors{f}, x);
  end
end



converged = 0;
iter = 1;
var_prod = cell(1, nvars);
fac_prod = cell(1, nfactors);

while ~converged & (iter <= engine.max_iter)
  if verbose, fprintf('iter %d\n', iter);  end
  
  % absorb
  old_var_prod = var_prod;
  for x=1:nvars
    var_prod{x} = mk_initial_pot(pot_type, x, ns, cnodes, onodes);
    for f=engine.fgraph.dep{x}
      var_prod{x} = multiply_by_pot(var_prod{x}, msg_fac_to_var{f,x});
    end
  end
  for f=1:nfactors
    fac_prod{f} = mk_initial_pot(pot_type, dom{f}, ns, cnodes, onodes);
    for x=dom{f}
      fac_prod{f} = multiply_by_pot(fac_prod{f}, msg_var_to_fac{x,f});
    end
  end

  % send msgs to neighbors
  old_msg_var_to_fac = msg_var_to_fac;
  old_msg_fac_to_var = msg_fac_to_var;
  converged = 1;
  for x=1:nvars
    %if verbose, disp(['var ' num2str(x) ' sending to fac ' num2str(engine.fgraph.dep{x})]); end
    for f=engine.fgraph.dep{x}
      temp = divide_by_pot(var_prod{x}, old_msg_fac_to_var{f,x});
      msg_var_to_fac{x,f} = normalize_pot(temp);
      if ~approxeq_pot(msg_var_to_fac{x,f}, old_msg_var_to_fac{x,f}, engine.tol), converged = 0; end
    end
  end
  for f=1:nfactors
    %if verbose, disp(['fac ' num2str(f) ' sending to var ' num2str(dom{f})]); end
    for x=dom{f}
      temp = divide_by_pot(fac_prod{f}, old_msg_var_to_fac{x,f});
      temp2 = multiply_by_pot(factors{f}, temp);
      temp3 = marginalize_pot(temp2, x, maximize);
      msg_fac_to_var{f,x} = normalize_pot(temp3);
      if ~approxeq_pot(msg_fac_to_var{f,x}, old_msg_fac_to_var{f,x}, engine.tol), converged = 0; end
    end
  end

  if iter==1
    converged = 0;
  end
  iter = iter + 1;
end

niter = iter - 1;
engine.niter = niter;

for x=1:nvars
  engine.marginal_nodes{x} = normalize_pot(var_prod{x});
end


