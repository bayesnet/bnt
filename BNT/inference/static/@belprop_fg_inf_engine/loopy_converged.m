function niter = loopy_converged(engine)
% LOOPY_CONVERGED Did loopy belief propagation converge? 0 means no, eles we return the num. iterations.
% function niter = loopy_converged(engine)
%
% We use a simple heuristic: we say convergence occurred if the number of iterations
% used was less than the maximum allowed.

if engine.niter == engine.max_iter
  niter = 0;
else
  niter = engine.niter;
end
