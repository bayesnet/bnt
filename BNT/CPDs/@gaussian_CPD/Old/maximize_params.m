function CPD = maximize_params(CPD, temp)
% MAXIMIZE_PARAMS Set the params of a CPD to their ML values (Gaussian)
% CPD = maximize_params(CPD, temperature)
%
% Temperature is currently only used for entropic prior on Sigma

% For details, see "Fitting a Conditional Gaussian Distribution", Kevin Murphy, tech. report,
% 1998, available at www.cs.berkeley.edu/~murphyk/papers.html
% Refering to table 2, we use equations 1/2 to estimate the covariance matrix in the untied/tied case,
% and equation 9 to estimate the weight matrix and mean.
% We do not implement spherical Gaussians - the code is already pretty complicated!

if ~adjustable_CPD(CPD), return; end

%assert(approxeq(CPD.nsamples, sum(CPD.Wsum)));
assert(~any(isnan(CPD.WXXsum)))
assert(~any(isnan(CPD.WXYsum)))
assert(~any(isnan(CPD.WYYsum)))

[self_size cpsize dpsize] = size(CPD.weights);

% Append 1s to the parents, and derive the corresponding cross products.
% This is used when estimate the means and weights simultaneosuly,
% and when estimatting Sigma.
% Let x2 = [x 1]'
XY = zeros(cpsize+1, self_size, dpsize); % XY(:,:,i) = sum_l w(l,i) x2(l) y(l)' 
XX = zeros(cpsize+1, cpsize+1, dpsize); % XX(:,:,i) = sum_l w(l,i) x2(l) x2(l)' 
YY = zeros(self_size, self_size, dpsize); % YY(:,:,i) = sum_l w(l,i) y(l) y(l)' 
for i=1:dpsize
  XY(:,:,i) = [CPD.WXYsum(:,:,i) % X*Y
	       CPD.WYsum(:,i)']; % 1*Y
  % [x  * [x' 1]  = [xx' x
  %  1]              x'  1]
  XX(:,:,i) = [CPD.WXXsum(:,:,i) CPD.WXsum(:,i);
	       CPD.WXsum(:,i)'   CPD.Wsum(i)];
  YY(:,:,i) = CPD.WYYsum(:,:,i);
end

w = CPD.Wsum(:);
% Set any zeros to one before dividing
% This is valid because w(i)=0 => WYsum(:,i)=0, etc
w = w + (w==0);

if CPD.clamped_mean
  % Estimating B2 and then setting the last column (the mean) to the clamped mean is *not* equivalent
  % to estimating B and then adding the clamped_mean to the last column.
  if ~CPD.clamped_weights
    B = zeros(self_size, cpsize, dpsize);
    for i=1:dpsize
      if det(CPD.WXXsum(:,:,i))==0
	B(:,:,i) = 0;
      else
	% Eqn 9 in table 2 of TR
	%B(:,:,i) = CPD.WXYsum(:,:,i)' * inv(CPD.WXXsum(:,:,i));
	B(:,:,i) = (CPD.WXXsum(:,:,i) \ CPD.WXYsum(:,:,i))';
      end
    end
    %CPD.weights = reshape(B, [self_size cpsize dpsize]);
    CPD.weights = B;
  end
elseif CPD.clamped_weights % KPM 1/25/02
  if ~CPD.clamped_mean % ML estimate is just sample mean of the residuals
    for i=1:dpsize
      CPD.mean(:,i) = (CPD.WYsum(:,i) - CPD.weights(:,:,i) * CPD.WXsum(:,i)) / w(i);
    end
  end
else % nothing is clamped, so estimate mean and weights simultaneously
  B2 = zeros(self_size, cpsize+1, dpsize);
  for i=1:dpsize
    if det(XX(:,:,i))==0  % fix by U. Sondhauss 6/27/99
      B2(:,:,i)=0;          
    else                    
      % Eqn 9 in table 2 of TR
      %B2(:,:,i) = XY(:,:,i)' * inv(XX(:,:,i));
      B2(:,:,i) = (XX(:,:,i) \ XY(:,:,i))';
    end                   
    CPD.mean(:,i) = B2(:,cpsize+1,i);
    CPD.weights(:,:,i) = B2(:,1:cpsize,i);
  end
end

% Let B2 = [W mu]
if cpsize>0
  B2(:,1:cpsize,:) = reshape(CPD.weights, [self_size cpsize dpsize]);
end
B2(:,cpsize+1,:) = reshape(CPD.mean, [self_size dpsize]);

% To avoid singular covariance matrices,
% we use the regularization method suggested in "A Quasi-Bayesian approach to estimating
% parameters for mixtures of normal distributions", Hamilton 91.
% If the ML estimate is Sigma = M/N, the MAP estimate is (M+gamma*I) / (N+gamma),
% where gamma >=0 is a smoothing parameter (equivalent sample size of I prior)

gamma = CPD.cov_prior_weight;

if ~CPD.clamped_cov
  if CPD.cov_prior_entropic % eqn 12 of Brand AI/Stat 99
    Z = 1-temp;
    % When temp > 1, Z is negative, so we are dividing by a smaller
    % number, ie. increasing the variance.
  else
    Z = 0;
  end
  if CPD.tied_cov
    S = zeros(self_size, self_size);
    % Eqn 2 from table 2 in TR
    for i=1:dpsize
      S = S + (YY(:,:,i) - B2(:,:,i)*XY(:,:,i));
    end
    %denom = max(1, CPD.nsamples + gamma + Z);
    denom = CPD.nsamples + gamma + Z;
    S = (S + gamma*eye(self_size)) / denom;
    if strcmp(CPD.cov_type, 'diag')
      S = diag(diag(S));
    end
    CPD.cov = repmat(S, [1 1 dpsize]);
  else 
    for i=1:dpsize      
      % Eqn 1 from table 2 in TR
      S = YY(:,:,i) - B2(:,:,i)*XY(:,:,i);
      %denom = max(1, w(i) + gamma + Z); % gives wrong answers on mhmm1
      denom = w(i) + gamma + Z;
      S = (S + gamma*eye(self_size)) / denom;
      CPD.cov(:,:,i) = S;
    end
    if strcmp(CPD.cov_type, 'diag')
      for i=1:dpsize      
	CPD.cov(:,:,i) = diag(diag(CPD.cov(:,:,i)));
      end
    end
  end
end


check_covars = 0;
min_covar = 1e-5;
if check_covars % prevent collapsing to a point
  for i=1:dpsize
    if min(svd(CPD.cov(:,:,i))) < min_covar
      disp(['resetting singular covariance for node ' num2str(CPD.self)]);
      CPD.cov(:,:,i) = CPD.init_cov(:,:,i);
    end
  end
end



