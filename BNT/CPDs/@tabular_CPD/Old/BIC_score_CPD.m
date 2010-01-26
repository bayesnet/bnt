function score = BIC_score_CPD(CPD, fam, data, ns, cnodes)
% BIC_score_CPD Compute the BIC score of a tabular CPD
% score = BIC_score_CPD(CPD, fam, data, ns, cnodes)

if iscell(data)
  local_data = cell2num(data(fam,:));
else
  local_data = data(fam, :);
end
counts = compute_counts(local_data, CPD.sizes);
CPT = mk_stochastic(counts); % MLE
tiny = exp(-700); 
CPT = CPT + (CPT==0)*tiny;  % replace 0s by tiny
LL = sum(log(CPT(:)) .* counts(:));
N = size(data, 2);
score = LL - 0.5*CPD.nparams*log(N);

