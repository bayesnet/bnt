% function [L,Ph,LL]=ffa(X,K,cyc,tol);
% 
% Fast Maximum Likelihood Factor Analysis using EM
%
% X - data matrix
% K - number of factors
% cyc - maximum number of cycles of EM (default 100)
% tol - termination tolerance (prop change in likelihood) (default 0.0001)
%
% L - factor loadings 
% Ph - diagonal uniquenesses matrix
% LL - log likelihood curve
%
% Iterates until a proportional change < tol in the log likelihood 
% or cyc steps of EM 
%

function [L,Ph,LL]=ffa(X,K,cyc,tol);

if nargin<4  tol=0.0001; end;
if nargin<3  cyc=100; end;

N=length(X(:,1));
D=length(X(1,:));
tiny=exp(-700);

X=X-ones(N,1)*mean(X);
XX=X'*X/N;
diagXX=diag(XX);

randn('seed', 0);
cX=cov(X);
scale=det(cX)^(1/D);
L=randn(D,K)*sqrt(scale/K);
Ph=diag(cX);

I=eye(K);

lik=0; LL=[];

const=-D/2*log(2*pi);


for i=1:cyc;

  %%%% E Step %%%%
  Phd=diag(1./Ph);
  LP=Phd*L;
  MM=Phd-LP*inv(I+L'*LP)*LP';
  dM=sqrt(det(MM));
  beta=L'*MM;
  XXbeta=XX*beta';
  EZZ=I-beta*L +beta*XXbeta;

  %%%% Compute log likelihood %%%%
  
  oldlik=lik;
  lik=N*const+N*log(dM)-0.5*N*sum(diag(MM*XX));
  fprintf('cycle %i lik %g \n',i,lik);
  LL=[LL lik];
  
  %%%% M Step %%%%

  L=XXbeta*inv(EZZ);
  Ph=diagXX-diag(L*XXbeta');

  if (i<=2)    
    likbase=lik;
  elseif (lik<oldlik)     
    disp('VIOLATION');
  elseif ((lik-likbase)<(1+tol)*(oldlik-likbase)||~isfinite(lik))  
    break;
  end;

end
