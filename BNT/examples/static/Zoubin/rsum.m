% row sum
% function Z=rsum(X)

function Z=rsum(X)

[N M]=size(X);

Z=zeros(N,1);

if M==1,
  Z=X;
elseif M<2*N,
  for m=1:M,
    Z=Z+X(:,m);
  end;
else
  for n=1:N
    Z(n)=sum(X(n,:));
  end;
end
