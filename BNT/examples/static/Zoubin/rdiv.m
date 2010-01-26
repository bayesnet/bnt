% function Z=rdiv(X,Y)
%
% row division: Z = X / Y row-wise
% Y must have one column 

function Z=rdiv(X,Y)

[N M]=size(X);
[K L]=size(Y);
if(N ~= K | L ~=1)
  disp('Error in RDIV');
  return;
end

Z=zeros(N,M);

if M<N,
  for m=1:M
    Z(:,m)=X(:,m)./Y;
  end
else
  for n=1:N
    Z(n,:)=X(n,:)/Y(n);
  end;
end;