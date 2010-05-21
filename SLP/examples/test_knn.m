function [ratio, ratiominus, ratioplus, tps]  = test_knn(BD,BDT,class,K)
% [epsilon, tps]  = test_knn(BD,BDT,class,K)
%
% epsilon = good classification rate on BDT,
% tps = computation time,
% BD = examples dataset using to find neighbors,
% BDT = set of examples to class,
% class = classification attribute number,
% K = vote on the K nearest neighbor
%
% francois.olivier.c.h@gmail.com

[N L]=size(BD);
Ltest=size(BDT,2);
fprintf(' test dataset size %d\r\n',Ltest);
E=setdiff(1:N,class);
place=0;
%for k=K
k=K;
 place=place+1;
 good=0;tic
 for i=1:Ltest
  %fprintf(' %d\r',i);
  res=knn(BD(E,:)',BD(class,:)',[1 2],BDT(E,i)',k);
  if res==BDT(class,i)
    good=good+1;
  end
 end
 tps(place)=toc;
 epsilon(1,place)=k;
 epsilon(2,place)=good/(Ltest);
%end

ratio = epsilon(2,place);
[ratiominus ratioplus] = confiance(ratio,L);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�
function [ypred]=knn(xapp,yapp,valY,X,k)
% knn implementation
%
% USE : [ypred]=knn(xapp,yapp,valY,X,k)
%
% Vincent Guigue 08/01/03

if nargin<4
  error('too few argumemnts');
elseif nargin<5
  k=3;
else
  if mod(k,2)==0
    error('k must be odd');
  end
end

if size(xapp,2)~=size(X,2)
  error('dimension incompatibility');
end
ndim = size(xapp,2);
nptxapp = size(xapp,1);
nptX = size(X,1);

% distance de X a xapp :
mat1 =  repmat(xapp, nptX,1);
%mat21 = reshape(X',1,nptX*ndim)
mat22 = repmat(X,1,nptxapp)';
mat2 = reshape(mat22 ,ndim, nptxapp*nptX)';
distance = mat1 - mat2 ;
distance = sum(distance.^2,2);
distance = reshape(distance,nptxapp,nptX);
[val kppv] = sort(distance,1);
% bilan sur les k premieres lignes
kppv = reshape(kppv(1:k,:),k*nptX,1);
Ykppv = yapp(kppv,1);
Ykppv = reshape(Ykppv,k,nptX);
% trouver le plus de reponses identique par colonne
vote = [];
for i=1:nptX
  for j=1:length(valY)
    vote(j,i)=size(find(Ykppv(:,i)==valY(j)),1);
  end
end
[val ind]=max(vote,[],1);
ypred = valY(ind);