function [t,nk] = minspan(IJC)
%MINSPAN Minimum weight spanning tree using Kruskal algorithm.
%[t,nk] = minspan(IJC)
%   IJC = n x 3 matrix arc list [i j c] of arc heads, tails, and costs
%     t = n-element logical vector, where 
%         t(i) = 1, if IJC(i,:) arc in spanning tree
%         t(i) = k, if IJC(i,:) arc in component k of forest
%    nk = number of components

% Copyright (c) 1998-2001 by Michael G. Kay
% Matlog Version 5 22-Aug-2001

% Input Error Checking ******************************************************
[n,cIJC] = size(IJC);
if cIJC ~= 3
   error('''IJC'' must be a 3-column matrix.')
elseif n < 1
   error('''IJC'' must have at least one row.')
elseif any(IJC(:,1) < 1) | any(any(~isint(IJC(:,[1 2]))))
   error('Invalid arc index in IJC.')
end
% End (Input Error Checking) ************************************************

i = IJC(:,1); j = abs(IJC(:,2));
m = max(max([i j]));

sidxIJ = argsort(IJC(:,3));
i = i(sidxIJ); j = j(sidxIJ);

t = logical(zeros(n,1));
k = 1;            % Current arc
nt = 0;           % Number of arcs in spanning tree
v = (1:m)';       % Arc labels

while nt < m - 1 & k <= n
   if (v(i(k)) ~= v(j(k)))
      v(v==v(j(k))) = v(i(k));
      t(k) = 1;
      nt = nt + 1;
   end
   k = k + 1;
end

idxIJ = invperm(sidxIJ);
t = t(idxIJ); i = i(idxIJ); j = j(idxIJ);

c = unique(v(unique([i; j])));   % Unique labels of arc vertices
nk = length(c);
if ~any(t), nk = 0; end          % Self-loop not a component

if nk > 1
   for k = 1:nk
      t(t~=0 & v(i)==c(k)) = k;  % Relabel to consecutive component numbers
   end
end

