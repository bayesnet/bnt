function [c,v] = scc(a,tol)

%       Finds the strongly connected sets of vertices
%                in the DI-rected G-raph of A
%          c = 0-1 matrix displaying accessibility
%          v = displays the equivalent classes
%
% v(i,j) is the j'th member of the i'th equiv class (0 padded)
%
% http://www.math.wsu.edu/math/faculty/tsat/matlab.html

[m,n] = size(a);
if m~=n 'Not a Square Matrix', break, end
b=abs(a); o=ones(size(a)); x=zeros(1,n);
msg='The Matrix is Irreducible !'; v='Connected Directed Graph !';
if (nargin==1) tol=n*eps*norm(a,'inf'); end

% Create a companion matrix
c = b>tol*o;
if (c==o)
  %  msg, break
  v = 1:length(a);
  return
end


% Compute accessibility in at most n-step paths
for k=1:n
  for j=1:n
    for i=1:n
      % If index i accesses j, where can you go ?
      if c(i,j) > 0  c(i,:) = c(i,:)+c(j,:); end
    end
  end
end
% Create a 0-1 matrix with the above information
c>zeros(size(a)); c=ans; if (c==o) msg, break, end

% Identify equivalence classes
d=c.*c'+eye(size(a)); d>zeros(size(a)); d=ans;
v=zeros(size(a));
for i=1:n find(d(i,:)); ans(n)=0; v(i,:)=ans; end

% Eliminate displaying of identical rows
i=1;
while(i<n)
  for k=i+1:n
    if v(k,1) == v(i,1)
      v(k,:)=x;
    end
  end
  i=i+1;
end
j=1;
for i=1:n
  if v(i,1)>0
    h(j,:)=v(i,:);
    j=j+1;
  end
end
v=h;



