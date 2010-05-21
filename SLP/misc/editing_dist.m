function D = editing_dist(dag1, dag2)
% D = editing_dist(dag1, dag2)
%
% d = 1 if arc1 <> arc2
%

[n1 m1]=size(dag1);
[n2 m2]=size(dag2);

if n1~=m1 | n1~=n2 | n2~=m2 | m1~=m2
  error('formats non compatibles');
  D=-inf;
end

de=abs(dag1-dag2);
mauv=find(triu(or(de,de')));
D=size(mauv,1);
