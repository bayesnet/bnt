function score = bayesian_score_CPD(CPD, local_ev)
% bayesian_score_CPD Compute the Bayesian score of a tabular CPD using uniform Dirichlet prior
% score = bayesian_score_CPD(CPD, local_ev)
%
% The Bayesian score is the log marginal likelihood

if iscell(local_ev)
 data = num2cell(local_ev);
else
 data =	local_ev;
end

score = dirichlet_score_family(compute_counts(data, CPD.sizes));
