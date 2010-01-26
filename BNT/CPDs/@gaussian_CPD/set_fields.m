function CPD = set_fields(CPD, varargin)
% SET_PARAMS Set the parameters (fields) for a gaussian_CPD object
% CPD = set_params(CPD, name/value pairs)
%
% The following optional arguments can be specified in the form of name/value pairs:
%
% mean       - mu(:,i) is the mean given Q=i
% cov        - Sigma(:,:,i) is the covariance given Q=i 
% weights    - W(:,:,i) is the regression matrix given Q=i 
% cov_type   - if 'diag', Sigma(:,:,i) is diagonal 
% tied_cov   - if 1, we constrain Sigma(:,:,i) to be the same for all i
% clamp_mean - if 1, we do not adjust mu(:,i) during learning 
% clamp_cov  - if 1, we do not adjust Sigma(:,:,i) during learning 
% clamp_weights - if 1, we do not adjust W(:,:,i) during learning
% clamp      - if 1, we do not adjust any params
% cov_prior_weight - weight given to I prior for estimating Sigma
% cov_prior_entropic - if 1, we also use an entropic prior for Sigma [0]
%
% e.g., CPD = set_params(CPD, 'mean', [0;0])

args = varargin;
nargs = length(args);
for i=1:2:nargs
  switch args{i},
   case 'mean',        CPD.mean = args{i+1}; 
   case 'cov',         CPD.cov = args{i+1}; 
   case 'weights',     CPD.weights = args{i+1}; 
   case 'cov_type',    CPD.cov_type = args{i+1}; 
   %case 'tied_cov',    CPD.tied_cov = strcmp(args{i+1}, 'yes');
   case 'tied_cov',    CPD.tied_cov = args{i+1};
   case 'clamp_mean',  CPD.clamped_mean = args{i+1};
   case 'clamp_cov',   CPD.clamped_cov = args{i+1};
   case 'clamp_weights',  CPD.clamped_weights = args{i+1};
   case 'clamp',  clamp = args{i+1};
    CPD.clamped_mean = clamp;
    CPD.clamped_cov = clamp;
    CPD.clamped_weights = clamp;
   case 'cov_prior_weight',  CPD.cov_prior_weight = args{i+1};
   case 'cov_prior_entropic',  CPD.cov_prior_entropic = args{i+1};
   otherwise,  
    error(['invalid argument name ' args{i}]);
  end
end
