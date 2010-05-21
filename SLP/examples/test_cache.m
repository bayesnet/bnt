% INITIALISATIONS
%%%%%%%%%%%%%%%%%

clear
%load asia5000

N=8;L=100;
ns=2*ones(1,N);
%data=asiab(:,1:5000);
m=5000;

bnet=mk_asia2_bnet;
data = cell(N,m);
for l = 1:m, data(:,l) = sample_bnet(bnet); end
data=cell2mat(data);asiab=data;
fprintf('Complete data have been created.');

names={ 'A' , 'S' , 'T' , 'L' , 'B' , 'O' , 'X' , 'D' };

scoring_fn='bic';
%scoring_fn='bayesian';

% greedy search sur asia sans le cache
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

res=zeros(2,10);i=0;

for L=100:100:500
i=i+1;L

tps=cputime;
D1=learn_struct_gs(data, ns, zeros(N,N),'scoring_fn',scoring_fn);
tps=cputime-tps
%subplot(1,4,2);draw_graph(D1,names,ones(1,N),xx,yy);

% score BIC (ce sont des chiffres moyen sur ma machine pour environ 5 test par catégories, et idem pour les autres res)
% 100 exemples : 24s
% 500 exemples : 33s
% 1500 exemples : 49s
% 5000 exemples : 112s
% score Bayesien
% 100 exemples : 30
% 1500 exemples : 39
% 5000 exemples : 102

res(1,i)=tps;

% greedy search sur asia avec le cache
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%cache=sparse(L,N+3);
cache=score_init_cache(N,L);

% QD le cache est vide
tps2=cputime;
[D2 best2]=learn_struct_gs(data, ns, zeros(N,N),'scoring_fn',scoring_fn,'cache',cache);
tps2=cputime-tps2

res(2,i)=tps2;
end

%subplot(1,4,3);draw_graph(D2,names,ones(1,N),xx,yy);
% L=100
% score BIC
% 100 exemples : 28s
% taille du cache : 4860 octets pour 94 entrées au lieu de 8272 pour une matrice de taille 94*11.
% 1500 exemples : 45s
% 5000 exemples : 77s
% score Bayesien
% 100 exemples : 28
% 1500 exemples : 35
% 5000 exemples : 78

% QD le cache existe
%tps3=cputime;
%[D3 best3]=learn_struct_gs(data, ns, zeros(N,N),'scoring_fn',scoring_fn,'cache', cache2);
%tps3=cputime-tps3
%subplot(1,4,4);draw_graph(D3,names,ones(1,N),xx,yy);
% L=100
% score BIC
% 100 exemples : 37s qd le cache existe.
% 1500 exemples : 69s 62
% 5000 exemples : 129s
% score Bayesien
% 100 exemples : 27   }
% 1500 exemples : 66  } -> phenomène du a une taille de cache trop petite
% 5000 exemples : 127 }    ca se stabilise à partir de L=200 voir conclusion



L=400;
cache=zeros(L,N+3);
tps0=zeros(2,20);i=0;
for taille=250:250:2000
i=i+1;taille
cache=score_init_cache(N,L);
data=asiab(:,1:taille);
tps=cputime;
D=learn_struct_gs(data, ns, zeros(N,N),'scoring_fn',scoring_fn);
tps0(1,i)=cputime-tps;
tps=cputime;
[D best]=learn_struct_gs(data, ns, zeros(N,N),'scoring_fn',scoring_fn,'cache',cache);
tps0(2,i)=cputime-tps
end

% taille BD :  100   200   300   400   500   600   700   800   900  1000  1100  1200  1300  1400  1500  1600  1700  1800  1900  2000
% ss cache  : 24.1  24.7  29.8  51.1  31.1  36.8  38.0  30.3  48.7  44.6  47.9  46.2  46.9  51.6  50.3  50.3  56.5  50.9  51.6  48.6
% avec cache: 22.9  23.0  27.1  46.6  27.6  32.5  32.5  26.0  41.3  36.6  39.3  37.5  37.0  40.3  37.4  39.0  41.9  37.5  40.0  34.4



%%%%%%%%%%%%%%
% CONCLUSION %
%%%%%%%%%%%%%%

% Si la base de données est petite, il vaut mieux ne pas utiliser de cache.
% Celui-ci deviens utile à partir de la taille 500 environ...
% Au final on arrive presque à gagner 25% du tps avec 1500 exemple et 40% avec 5000 exemples sur le temps d'execution sur cet exemple.
% Utiliser un cache vide, ou un cache possédant deja des valeurs ne change pas grand chose, à part si celui-ci est de taille trop petite ^^(voir res pour 5000 exemples et score bic et res2 pour 1500 ex et score bayésien et res3 pour 500 ex et bic et res4 pour 1500 et bic)).
% Sinon qd le cache est trop gd cela n'est pas trop grave.
% En fait, utiliser une matrice classique prend moins de place mémoire qd le cache est plein.
%
% res= (L=100    200    300    400    500    600    700    800    900   1000)
% ssc   112
% cvide 86.5   63.6   64.4   65.1   63.9   64.1   64.4   64.5   64.2   64.8
% cplein127.7   63.4   63.8   63.7   63.4   63.9   63.5   63.9   63.6   64.2
%
% res2= (L=20     40     60     80    100    120    140    160    180    200    300    400    500    600)
% ssc    39
% cvide       46.5   49.2   49.7   52.3   41.9   36.3   46.5   47.1   47.3   46.2   46.4   34.9   35.4   35.9
% cplein 44.9   46.2   47.7   51.1   69.0   36.0   37.1   37.6   48.2   47.7   48.2   47.7   36.7   36.2
%
% res3= (L=50    100    150    200    250    300    350    400    450    500    550    600)
% ssc    33
% cvide  41.9   28.8   28.5   28.6   28.4   28.7   28.4   28.9   28.7   28.8   28.5   28.6
% cplein 39.1   46.3   28.9   29.3   29.0   29.5   29.1   29.6   29.3   29.0   29.4   29.3
%
% res4= (L=50    100    150    200    250    300    350    400)
% ssc    49
% cvide  55.3    39.0   37.2   39.2   37.9   38.3   37.8   38.4
%
%
% ON voit donc qu'une fois la taille du cache bien choisie, celui-ci ne fait pas perdre de perfs sur les petites BD,
% IL ne reste donc plus qu'a trouver une formule qui donne une taille de cache optimum en fonction du nombre de noeud,
% et esperer que celle-ci ne soit pas trop exponnentielle...




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ANNEXE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ChronoMetrage pour 1500 exemples.
% (pas à jour, c'était comme ca au debut maintenant le prétraitement est passé à ~~0.25s)
%%%%%%%%%%%%%%%

%                n'existe pas ds le cache  |  existe     |  ss cache
% -----------------------------------------|-------------|--------------
% Prétraitement   ==> 0.0044~64 secondes   | 0.0063~84 s |              <== trop long
% Lecture         ==> 0.0004~09 secondes   | 0.0004~09 s |
% if ~exist                                |             |
%   calcul        ==> 0.0076~87 secondes   |             | 0.0076~87 s  <== par rapport a ca
%   ecriture      ==> 0.0005~05 secondes   |             |
% _________________________________________|_____________|_______________
% TOTAL           === 1,3~1,6 ms           |  0,7~,9 ms  | 0,9 ms
%                        |
%                        |__> à évoluer vers ????

% enlever les commentaires sur les tic,toc en dur ds les fonctions

%  tic
%  [scorr cac]=score_family(2,[3 4 5],'tabular', 'bic', 2*ones(1,N),1:N, data, {}, cache);
%  toc
%  tic
%  [scorr cac]=score_family(2,[3 4 5],'tabular', 'bic', 2*ones(1,N),1:N, data, {}, cac);
%  toc
%  tic
%  scorr=score_family(2,[3 4 5],'tabular', 'bic', 2*ones(1,N),1:N, data, {});
%  toc
