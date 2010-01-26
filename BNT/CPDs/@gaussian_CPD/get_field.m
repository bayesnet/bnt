function val = get_params(CPD, name)
% GET_PARAMS Get the parameters (fields) for a gaussian_CPD object
% val = get_params(CPD, name)
%
% The following fields can be accessed
%
% mean       - mu(:,i) is the mean given Q=i
% cov        - Sigma(:,:,i) is the covariance given Q=i 
% weights    - W(:,:,i) is the regression matrix given Q=i 
%
% e.g., mean = get_params(CPD, 'mean')

switch name
 case 'mean',      val = CPD.mean;
 case 'cov',       val = CPD.cov;
 case 'weights',   val = CPD.weights;
 otherwise,
  error(['invalid argument name ' name]);
end
