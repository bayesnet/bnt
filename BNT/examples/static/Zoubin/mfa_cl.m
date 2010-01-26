% function [lik, likv]=mfa_cl(X,Lh,Ph,Mu,Pi);
% 
% Calculates log likelihoods of a data set under a mixture of factor
% analysis model.
%
% X - data matrix
% Lh - factor loadings 
% Ph - diagonal uniquenesses matrix
% Mu - mean vectors
% Pi - priors
%
% lik - log likelihood of X 
% likv - vector of log likelihoods
% 
% If 0 or 1 output arguments requested, lik is returned. If 2 output
% arguments requested, [lik likv] is returned.

function [lik, likv]=mfa_cl(X,Lh,Ph,Mu,Pi);

N=length(X(:,1));
D=length(X(1,:));
K=length(Lh(1,:));
M=length(Pi);

if (abs(sum(Pi)-1) > 1e-6) 
  disp('ERROR: Pi should sum to 1');
  return;
elseif ((size(Lh) ~= [D*M K]) | (size(Ph) ~= [D 1]) | (size(Mu) ~= [M D]) ...
  | (size(Pi) ~= [M 1] & size(Pi) ~= [1 M]))   
  disp('ERROR in input matrix sizes');
  return;
end;  

tiny=exp(-744);
const=(2*pi)^(-D/2);

I=eye(K);
Phi=1./Ph;
Phid=diag(Phi);
for k=1:M  
  Lht=Lh((k-1)*D+1:k*D,:);
  LP=Phid*Lht;
  MM=Phid-LP*inv(I+Lht'*LP)*LP';
  dM=sqrt(det(MM));      	
  Xk=(X-ones(N,1)*Mu(k,:)); 
  XM=Xk*MM; 
  H(:,k)=const*Pi(k)*dM*exp(-0.5*sum((XM.*Xk)'))'; 	
end;

Hsum=rsum(H); 				

likv=log(Hsum+(Hsum==0)*tiny);
lik=sum(likv);

