% function [Lh,Ph,Mu,Pi,LL]=mfa(X,M,K,cyc,tol);
% 
% Maximum Likelihood Mixture of Factor Analysis using EM
%
% X - data matrix
% M - number of mixtures (default 1)
% K - number of factors in each mixture (default 2)
% cyc - maximum number of cycles of EM (default 100)
% tol - termination tolerance (prop change in likelihood) (default 0.0001)
%
% Lh - factor loadings 
% Ph - diagonal uniquenesses matrix
% Mu - mean vectors
% Pi - priors
% LL - log likelihood curve
%
% Iterates until a proportional change < tol in the log likelihood 
% or cyc steps of EM 

function [Lh, Ph,  Mu, Pi, LL] = mfa(X,M,K,cyc,tol)

if nargin<5   tol=0.0001; end;
if nargin<4   cyc=100; end;
if nargin<3   K=2; end;
if nargin<2   M=1; end;

N=length(X(:,1));
D=length(X(1,:));
tiny=exp(-700);

%rand('state',0);

fprintf('\n');

if (M==1)
  [Lh,Ph,LL]=ffa(X,K,cyc,tol);
  Mu=mean(X);
  Pi=1;
else
  if N==1
    mX = X;
  else
    mX=mean(X);
  end
  cX=cov(X);
  scale=det(cX)^(1/D);
  randn('state',0); 
  Lh=randn(D*M,K)*sqrt(scale/K);
  Ph=diag(cX)+tiny;
  Pi=ones(M,1)/M;
  %randn('state',0); 
  Mu=randn(M,D)*sqrtm(cX)+ones(M,1)*mX;
  oldMu=Mu;
  I=eye(K);

  lik=0;
  LL=[];

  H=zeros(N,M); 	% E(w|x) 
  EZ=zeros(N*M,K);
  EZZ=zeros(K*M,K);
  XX=zeros(D*M,D);
  s=zeros(M,1);
  const=(2*pi)^(-D/2);
  %%%%%%%%%%%%%%%%%%%%
  for i=1:cyc;

    %%%% E Step %%%%

    Phi=1./Ph;
    Phid=diag(Phi);
    for k=1:M
      Lht=Lh((k-1)*D+1:k*D,:);
      LP=Phid*Lht;
      MM=Phid-LP*inv(I+Lht'*LP)*LP';
      dM=sqrt(det(MM));      	
      Xk=(X-ones(N,1)*Mu(k,:)); 
      XM=Xk*MM;
      H(:,k)=const*Pi(k)*dM*exp(-0.5*rsum(XM.*Xk)); 	
      EZ((k-1)*N+1:k*N,:)=XM*Lht;
    end;
    
    Hsum=rsum(H);
    oldlik=lik;
    lik=sum(log(Hsum+(Hsum==0)*exp(-744)));

    Hzero=(Hsum==0); Nz=sum(Hzero); 
    H(Hzero,:)=tiny*ones(Nz,M)/M; 
    Hsum(Hzero)=tiny*ones(Nz,1);
    
    H=rdiv(H,Hsum); 				
    s=csum(H);
    s=s+(s==0)*tiny;
    s2=sum(s)+tiny;
    
    for k=1:M  
      kD=(k-1)*D+1:k*D;
      Lht=Lh(kD,:);
      LP=Phid*Lht;
      MM=Phid-LP*inv(I+Lht'*LP)*LP';
      Xk=(X-ones(N,1)*Mu(k,:)); 
      XX(kD,:)=rprod(Xk,H(:,k))'*Xk/s(k); 
      beta=Lht'*MM;
      EZZ((k-1)*K+1:k*K,:)=I-beta*Lht +beta*XX(kD,:)*beta'; 
    end;

    %%%% log likelihood %%%%

    LL=[LL lik];
    fprintf('cycle %g   \tlog likelihood %g ',i,lik);
    
    if (i<=2)
      likbase=lik;
    elseif (lik<oldlik) 
      fprintf(' violation');
    elseif ((lik-likbase)<(1 + tol)*(oldlik-likbase)||~isfinite(lik)) 
      break;
    end;

    fprintf('\n');
    
    %%%% M Step %%%%
    
    % means and covariance structure
    
    Ph=zeros(D,1);
    for k=1:M
      kD=(k-1)*D+1:k*D;
      kK=(k-1)*K+1:k*K;
      kN=(k-1)*N+1:k*N;

      T0=rprod(X,H(:,k));
      T1=T0'*[EZ(kN,:) ones(N,1)];
      XH=EZ(kN,:)'*H(:,k);
      T2=inv([s(k)*EZZ(kK,:) XH; XH' s(k)]);
      T3=T1*T2;
      Lh(kD,:)=T3(:,1:K);
      Mu(k,:)=T3(:,K+1)';
      T4=diag(T0'*X-T3*T1')/s2;
      Ph=Ph+T4.*(T4>0); 
    end;

    Phmin=exp(-700);
    Ph=Ph.*(Ph>Phmin)+(Ph<=Phmin)*Phmin; % to avoid zero variances

    % priors
    Pi=s'/s2;
    
  end;
  fprintf('\n');
end;


