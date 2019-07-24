function score = cond_mutual_info_score(i,si,j,sj,c,sc,data)
% G = cond_mutual_info_score(i,si,j,sj,c,sc,data)
% Only for tabular node which values are 1,2,...,size .
% si is size of node i, sj is size of node j, sc is the size of node c.
% data(i,m) is the node i in the case m.
% 
%
% pphilippe.leray@univ-nantes.fr, francois.olivier.c.h@gmail.com

[n N]=size(data);
Pc=hist(data(c,:),1:sc)/N;
score=0;

for cvalue=1:sc,
    ind=find(data(c,:)==cvalue);
    Nj=hist(data(j,ind),1:sj);
    Ni=hist(data(i,ind),1:si);
    NiNj=Ni'*Nj;

    for k=1:si
        ind2=find(data(i,ind)==k) ;
        Nij(k,:) = hist(data(j,ind(ind2)),1:sj);
    end

    % sommons les valeurs non-infinies:
    ind=find(NiNj~=0 & Nij~=0);
    score=score+Pc(cvalue)*sum(sum(Nij(ind).*log(N*Nij(ind)./NiNj(ind))/N));
end
