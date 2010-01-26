function [engine, ll, niter] = enter_soft_evidence(engine, local_evidence)
% ENTER_SOFT_EVIDENCE Propagate evidence using belief propagation
% [engine, ll, niter] = enter_soft_evidence(engine, local_evidence)
%
% local_evidence{i}(j) = Pr(observation at node i | S(i)=j)
%
% The log-likelihood is not computed; ll = 0.
% niter contains the number of iterations used 

ll = 0;
mrf2 = engine.mrf2;
[bel, niter] = bp_mrf2(mrf2.adj_mat, mrf2.pot, local_evidence, ...
		       'max_iter', engine.max_iter, 'momentum', engine.momentum, ...
		       'tol', engine.tol, 'maximize', 0, 'verbose', engine.verbose);
engine.bel = bel;
