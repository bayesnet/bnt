function mpe = find_mpe(engine, local_evidence)
% FIND_MPE Find the most probable explanation of the data  
% function mpe = find_mpe(engine, local_evidence
%
% local_evidence{i}(j) = Pr(observation at node i | S(i)=j)
%
% This finds the marginally most likely value for each hidden node.
% It may give inconsistent results if there are ties.

[mpe, niter] = bp_mpe_mrf2(engine.mrf2.adj_mat, engine.mrf2.pot, local_evidence, ...
			   'max_iter', engine.max_iter, 'momentum', engine.momentum, ...
			   'tol', engine.tol);
