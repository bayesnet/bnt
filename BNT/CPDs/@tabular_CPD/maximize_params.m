function CPD = maximize_params(CPD, temp)
% MAXIMIZE_PARAMS Set the params of a tabular node to their ML/MAP values.
% CPD = maximize_params(CPD, temp)

if ~adjustable_CPD(CPD), return; end

%assert(approxeq(sum(CPD.counts(:)), CPD.nsamples)); % false!
switch CPD.prior_type
 case 'none',
  counts = reshape(CPD.counts, size(CPD.CPT));
  CPD.CPT = mk_stochastic(counts);
 case 'dirichlet',
  counts = reshape(CPD.counts, size(CPD.CPT));
  CPD.CPT = mk_stochastic(counts + CPD.dirichlet);
 
 % case 'entropic',
%   % For an HMM,
%   % CPT(i,j) = pr(X(t)=j | X(t-1)=i) = transprob(i,j)
%   % counts(i,j) = E #(X(t-1)=i, X(t)=j) = exp_num_trans(i,j)
%   Z = 1-temp;
%   fam_sz = CPD.sizes;
%   psz = prod(fam_sz(1:end-1));
%   ssz = fam_sz(end);
%   counts = reshape(CPD.counts, psz, ssz);
%   CPT = zeros(psz, ssz);
%   for i=CPD.entropic_pcases(:)'
%     [CPT(i,:), logpost] = entropic_map_estimate(counts(i,:), Z);
%   end
%   non_entropic_pcases = mysetdiff(1:psz, CPD.entropic_pcases);
%   for i=non_entropic_pcases(:)'
%     CPT(i,:) = mk_stochastic(counts(i,:));
%   end
%   %for i=1:psz
%   %  [CPT(i,:), logpost] = entropic_map(counts(i,:), Z);
%   %end
%   if CPD.trim & (temp < 2) % at high temps, we would trim everything!
%     % grad(j) = d log lik / d theta(i ->j)
%     % CPT(i,j) = 0 => counts(i,j) = 0
%     % so we can safely replace 0s by 1s in the denominator
%     denom = CPT(i,:) + (CPT(i,:)==0);
%     grad = counts(i,:) ./ denom;
%     trim = find(CPT(i,:) <= exp(-(1/Z)*grad)); % eqn 32
%     if ~isempty(trim)
%       CPT(i,trim) = 0;
%       if all(CPD.trimmed_trans(i,trim)==0) % trimming for 1st time
% 	disp(['trimming CPT(' num2str(i) ',' num2str(trim) ')']) 
%       end
%       CPD.trimmed_trans(i,trim) = 1;
%     end
%   end
%   CPD.CPT = myreshape(CPT, CPD.sizes);
end
