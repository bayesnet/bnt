function [n,edges,nbedges,xechan] = hist_ic(x,crit)

%HIST_IC  optimal Histogram based on IC information criterion
%
%   [N,EDGES,NBEDGES,XECHAN] = HIST-IC(X,CRIT) 
%	bins the elements of X into an optimal number of bins according
%	to a cost function based on Akaike's Criterion.
%
%
%   CRIT = 1 | 2 | 3  (choose one of the 3 possible criterium)  (default 3)
%          4          (returns the initial histogram instead of the optimal one) 
%
%
%   N = cell array containing the distribution of each column of X
%   (or a vector if X is a column vector)
%   EDGES = cell array containing the bin edges of each column of X
%   (or a vector if X is a column vector)
%   NBEDGES = vector containing the number of bin edges for each column of X
%   (or a number if X is a column vector)
%   XECHAN = discretized version of X
%
%   Ref : O. Colot et al., Information Criteria and Abrupt Changes in
%         Probability Laws, Signal Processing VII: Theory and Applications
%	  pp.1855-18858, September 1994
%
%   F. El-Matouat, O. Colot 2000 (first version)
%   Revised 01-06-2001 by Ph. Leray - philippe.leray@univ-nantes.fr
%
%
%   Things to do :
%	* Call criteron by a name ('aic','xxx', ...) instead of a number
%


if nargin == 0
    error('Requires one or two arguments.')
end

if nargin == 1
    crit = 3;
end;

if min(size(x))==1, x = x(:); end

if isstr(x)
    error('Input argument must be numeric.')
end

if isempty(x),
	error('No elements to count')
end


[nb_l,nb_c]=size(x);

% Outputs declaration
xechan=zeros(nb_l,nb_c);

edges=cell(nb_c,1);
% Local variables
maxi = max(x);
mini = min(x);

%% Erreur ? ancien code :
%% nb_clas_ini=2*round(sqrt(nb_l)-1);	% article Fatima

nb_clas_ini=round(2*sqrt(nb_l)-1);	

pas_ini=(maxi-mini)/nb_clas_ini;	% initial step

for j=1:nb_c,

	% optimal histogram for each column of X
	histo_ini =hist(x(:,j),nb_clas_ini);		% initial histogram

	if (crit~=4)
		[hist_opt,pas_opt]=hist1_ic(histo_ini,nb_l,pas_ini(j),nb_clas_ini,crit);
	else
		fprintf('Histo initial\n');
		hist_opt=histo_ini;
		pas_opt=ones(1,nb_clas_ini)*pas_ini(j);
	end;
	nbedges(j)=size(hist_opt,2);
	edges{j}=mini(j)+cumsum(pas_opt(1:nbedges(j)-1)); %+1e-7;
        [n{j} xechan(:,j)]=histc(x(:,j),[-inf edges{j} inf]);
	n{j}=n{j}(1:end-1);
end

if (nb_c==1)
	n=n{1}; edges=edges{1};
end


% ============================== subfunctions

function [hist_opt,step_opt]=hist1_ic(histo,nb,step_ini,m,critere);

%HIST1_IC  optimal Histogram based on IC information criterion
%
%   [HIST_OPT, STEP_OPT] = HIST1_IC(HISTO, NB, STEP_INI, NBSTEP_INI, CRIT) 
%       fusion of an 1D histogramme (HISTO) according to an IC criterion (CRIT)
%
%   This function is mainly an internal function used by HIST_IC
%
%   Ref : O. Colot et al., Information Criteria and Abrupt Changes in
%         Probability Laws, Signal Processing VII: Theory and Applications
%         pp.1855-18858, September 1994
%
%   F. El-Matouat, O. Colot 2000 (first version)
%   Revised 11-06-2001 by Ph. Leray
%
%
%   Things to do :
%       * Call criteron by a name ('aic','xxx', ...) instead of a number
%


aic=[];
aic2=[];	

% Initialisation
histt=histo;
teta = histt/nb;
pas = step_ini*ones(1,m);

% Calcul de l'ensemble des histogrammes optimaux

for z=1:m

	% Calcul de AIC pour l'union entre hist(indice,u) et hist(indice,u+1)
	aic2 = [aic2 cal_aic(nb,teta,pas,m+1-z,critere)];

	if (z~=m)
		% Calcul des couples de classes adjacentes 
		if critere==1			
			penalite=(2*(m-z)-1)/nb;
		elseif critere==2
			penalite=(m-z-1)*(1+log(nb))/nb;
		else
			penalite=(m-z)*(1+log(log(nb)))/nb;
		end
 		aic=cla_adj(nb,teta,pas,step_ini,m-z+1,histt,penalite,aic);
		% Recherche de la valeur min du crit�re pour les classes adjacentes
		[min_aic classe]=min(aic(1:(m-z)));

	
		% Fusion de hist(classe) et hist(classe+1)
		nb_pas1=pas(classe)/step_ini;
		nb_pas2=pas(classe+1)/step_ini;

		ess=round( nb_pas1*histt(classe)+nb_pas2*histt(classe+1) );
		teta(classe)=ess / nb;			
		histt(classe)=ess / (nb_pas1+nb_pas2);
		pas(classe)=pas(classe)+pas(classe+1);


		% Cr�ation du nouvel histogramme
		itemp = setdiff(1:m+1-z,classe+1);
		histt = histt(itemp);
		pas = pas(itemp);
		teta = teta(itemp);
	end
end

% Recherche du crit�re minimun AIC
[min_AIC fusion]=min(aic2(1:m));

% Initialisation de histo
histt=histo;
teta = histt/nb;
pas = step_ini*ones(1,m);

% Calcul de l'histogramme optimal
		
for z=1:fusion-1

	% Calcul des couples de classes adjacentes 
	if critere==1			
		penalite=(2*(m-1)-1)/nb;
	elseif critere==2
		penalite=(m-2)*(1+log(nb))/nb;
	else
		penalite=(m-1)*(1+log(log(nb)))/nb;
	end

	aic=cla_adj(nb,teta,pas,step_ini,m,histt,penalite,aic);


	% Recherche de la valeur min du crit�re pour les classes adjacentes
	[min_aic classe]=min(aic(1:m-1));
			
	% Fusion de hist(indice,classe) et hist(indice,classe+1)

	nb_pas1=pas(classe)/step_ini;
	nb_pas2=pas(classe+1)/step_ini;
	
	teta(classe)=(round(nb_pas1*histt(classe)+nb_pas2*histt(classe+1)))/nb;
	histt(classe)=(nb_pas1*histt(classe)+nb_pas2*histt(classe+1))/(nb_pas1+nb_pas2);
	pas(classe)=pas(classe)+pas(classe+1);

	% Cr�ation du nouvel histogramme
				
	itemp=setdiff(1:m,classe+1);
	histt = histt(itemp);
	pas = pas(itemp);
	teta = teta(itemp);
	%aic=zeros(1,m-1);
				
	m=m-1;
end
hist_opt=histt;
step_opt=pas;


%=====================================================
% Calcul du Critere pour l'ensemble des classes

function akaike=cal_aic(size_ech,teta,pas,m,critere);


if critere==1
	a=(2*m-1)/size_ech;
elseif critere==2
	a=(m-1)*(1+log(size_ech))/size_ech;
else
	a=m*(1+log(log(size_ech)))/size_ech;
end

indu = find(teta);
akaike = a - 2*sum(teta(indu).*log(teta(indu)./pas(indu)));


%=====================================================
% Cla_adj.m
% aic=cla_adj(taille,indice,teta,pas,pas_ini,m,hist,penalite,aic)
% taille=nombre d'�l�ments dans chacune des hypotheses; 
% indice=numero de la classe;
% Calcul du critere de Akaike pour l'histogramme totale avec 
% fusion de deux classes adjacentes u et (u+1).

function aic=cla_adj(size_ech,teta,pas,pas_ini,m,hist,penalite,aic);


for u=1:m-1

	cumul=0;

	% This loop is faster than a sum of a vectorised computation !
	for x=1:m			
		if x~=u & x~=u+1 & teta(x)~=0
			cumul=cumul+teta(x)*log(teta(x)/pas(x));
		end			
	end
				
								
	nb_pas1=pas(u)/pas_ini;
	nb_pas2=pas(u+1)/pas_ini;

	b=( round(nb_pas1*hist(u)+nb_pas2*hist(u+1) ) ) / size_ech;

	if b~=0
		c=2*b*log( b / ( pas(u) + pas(u+1) ) );
	else
		c=0;
	end

	aic(u)=penalite-2*cumul-c;				
end
