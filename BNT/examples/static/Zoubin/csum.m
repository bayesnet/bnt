% column sum
% function Z=csum(X)

function Z=csum(X)

N=length(X(:,1));
if (N>1)
  Z=sum(X);
else
  Z=X;
end;