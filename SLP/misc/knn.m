function [ypred]=knn(xapp,yapp,valY,X,k)

%
% knn implementation
%
% USE : [ypred]=knn(xapp,yapp,valY,X,k)
%
% Vincent Guigue 08/01/03

% check nargin

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