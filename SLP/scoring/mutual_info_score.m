function score = mutual_info_score(i,si,j,sj,data)
% G = mutual_info_score(i,si,j,sj,data)
% Only for tabular node which values are 1,2,...,size .
% si is size of node i, sj is size of node j.
% data(i,m) is the node i in the case m.
% 
% Ref :
% C. Chow and C. Liu (1968). Approximating discrete probability distributions with dependence trees. 
% IEEE Transactions on Information Theory, 14(3):462--467, May 1968.
%
% francois.olivier.c.h@gmail.com, philippe.leray@univ-nantes.fr, wangxiangyang@sjtu.edu.cn

[n N]=size(data);
Nj=hist(data(j,:),1:sj);
Ni=hist(data(i,:),1:si);
NiNj=Ni'*Nj;

for k=1:si
 ind=find(data(i,:)==k) ;
 Nij(k,:) = hist(data(j,ind),1:sj);
end

% sommons les valeurs non-infinies:
ind=find(NiNj~=0 & Nij~=0);
score=sum(sum(Nij(ind).*log(N*Nij(ind)./NiNj(ind))/N));
