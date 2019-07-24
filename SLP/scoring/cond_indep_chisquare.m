function [CI, Chi2, alpha2] = cond_indep_chisquare(X, Y, S, Data, test, alpha, ns)
% COND_INDEP_CHISQUARE Test if X indep Y given Z
%                      using either chisquare test or likelihood ratio test G2
%
% [CI Chi2 Prob_Chi2] = cond_indep_chisquare(X, Y, S, Data, test, alpha, node_sizes)
%
% Input :
%       Data is the data matrix, N columns * NbVar rows
%       X is the index of variable X in Data matrix
%       Y is the index of variable Y in Data matrix
%       S are the indexes of variables in set S
%       alpha is the significance level (default: 0.05)
%       test = 'pearson' for Pearson's chi2 test
%		   'LRT' for G2 likelihood ratio test (default)
%       node_sizes (default: max(Data'))
%
% Output :
%       CI = test result (1=conditional independency, 0=no)
%       Chi2 = chi2 value (-1 if not enough data to perform the test --> CI=0)
%
%
% V1.4 : 24 july 2003 (Ph. Leray - philippe.leray@univ-nantes.fr)
%
%
% Things to do :
% - do not use 'find' in nij computation (when S=empty set)
% - find a better way than 'warning off/on' in tmpij, tmpijk computation
%

if nargin < 5, test = 'LRT'; end
if nargin < 6, alpha = 0.05; end
if nargin < 7, ns = max(Data'); end

Data=Data';

N = size(Data,1);
qi=ns(S);
tmp=[1 cumprod(qi(1:end-1))];
qs=1+(qi-1)*tmp';
if isempty(qs),
    nij=zeros(ns(X),ns(Y));
    df=prod(ns([X Y])-1)*prod(ns(S));
else

%   Commented by Mingyi
%    nijk=zeros(ns(X),ns(Y),qs);
%    tijk=zeros(ns(X),ns(Y),qs);
%   Commention ends
%   Added by Mingyi
    nijk=zeros(ns(X),ns(Y),1);
    tijk=zeros(ns(X),ns(Y),1);
%   Addition ends
    df=prod(ns([X Y])-1)*qs;
end


if (N<10*df)
    % Not enough data to perform the test
    Chi2=-1;
    CI=0;

elseif isempty(S)
    for i=1:ns(X),
        for j=1:ns(Y),
            nij(i,j)=length(find((Data(:,X)==i)&(Data(:,Y)==j))) ;
        end
    end
    restr=find(sum(nij,1)==0);
    if ~isempty(restr)
        nij=nij(:,find(sum(nij,1)));
    end

    tij=sum(nij,2)*sum(nij,1)/N ;

 switch test
    case 'pearson',
        tmpij=nij-tij;

        [xi yj]=find(tij<10);
        for i=1:length(xi),
           tmpij(xi(i),yj(i))=abs(tmpij(xi(i),yj(i)))-0.5;
        end

        warning off;
        tmp=(tmpij.^2)./tij;
        warning on;
        tmp(find(tmp==Inf))=0;

    case 'LRT',
        warning off;
        tmp=nij./tij;
        warning on;
        tmp(find(tmp==Inf | tmp==0))=1;
        tmp(find(tmp~=tmp))=1;
        tmp=2*nij.*log(tmp);

    otherwise,
        error(['unrecognized test ' test]);
    end

    Chi2=sum(sum(tmp));
    alpha2=1-chisquared_prob(Chi2,df);
    CI=(alpha2>=alpha) ;

else
    SizeofSSi=1;
    for exemple=1:N,
        i=Data(exemple,X);
        j=Data(exemple,Y);
        Si=Data(exemple,S)-1;
        %Added by Mingyi
        if exemple==1
            SSi(SizeofSSi,:)=Si;
            nijk(i,j,SizeofSSi)=1;
        else
            flag=0;
            for iii=1:SizeofSSi
                if isequal(SSi(iii,:),Si)
                    nijk(i,j,iii)=nijk(i,j,iii)+1;
                    flag=1;
                end
            end
            if flag==0
                SizeofSSi=SizeofSSi+1;
                SSi(SizeofSSi,:)=Si;
                nijk(i,j,SizeofSSi)=1;
            end
        end
        %Addition ends
        %Commented by Mingyi
%         k=1+Si*tmp';
%         nijk(i,j,k)=nijk(i,j,k)+1;
        %Commention ends
    end

    nik=sum(nijk,2);
    njk=sum(nijk,1);
    N2=sum(njk);

 %   for k=1:qs,         %Commented by Mingyi
    for k=1:SizeofSSi    %Added by Mingyi
        if N2(:,:,k)==0
            tijk(:,:,k)=0;
        else
            tijk(:,:,k)=nik(:,:,k)*njk(:,:,k)/N2(:,:,k);
        end
    end

    switch test
    case 'pearson',
        tmpijk=nijk-tijk;

        [xi yj]=find(tijk<10);
        for i=1:length(xi),
            tmpijk(xi(i),yj(i))=abs(tmpijk(xi(i),yj(i)))-0.5;
        end

        warning off;
        tmp=(tmpijk.^2)./tijk;
        warning on;
        tmp(find(tmp==Inf))=0;

    case 'LRT',
        warning off;
        tmp=nijk./tijk;
        warning on;
        tmp(find(tmp==Inf | tmp==0))=1;
        tmp(find(tmp~=tmp))=1;
        tmp=2*nijk.*log(tmp);

    otherwise,
        error(['unrecognized test ' test]);
    end

    Chi2=sum(sum(sum(tmp)));
    alpha2=1-chisquared_prob(Chi2,df);
    CI=(alpha2>=alpha) ;

end
clear tijk
clear nijk
clear nij
clear tij
clear tmpijk