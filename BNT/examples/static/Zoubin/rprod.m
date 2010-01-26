% row product
% function Z=rprod(X,Y)

function Z=rprod(X,Y)

if(length(X(:,1)) ~= length(Y(:,1)) | length(Y(1,:)) ~=1)
  disp('Error in RPROD');
  return;
end

Z=zeros(size(X));

for i=1:length(X(1,:))
  Z(:,i)=X(:,i).*Y;
end
