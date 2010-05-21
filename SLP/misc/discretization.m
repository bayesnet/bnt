function [appD, testD, bornes] = discretization(critere, continues, miss, app, test)
% [appD, testD, bornes] = discretization(critere, continious, miss, app, test)
%
% Inputs :
%    critere = 1, 2, 3 or 4 (see hict_ic for details)
%    continious = vector of continious variables to discretize
%    miss ~= 0 if it exists missing values coded by 'miss' value
%    app = Learning base
%    test = Test base (only the learning base is used to make the discretization rules) [optionnal]
%
% Outputs :
%    appD = Learning base with discretized entries on 'continious' variables
%    testD = Test base with discretized entries on 'continious' variables
%    bornes = limits of discretization intervals found by hist_ic
%

tt=1;
if nargin<5, test=[]; tt=0; end
app = app';
test = test';
[ma, Na] = size(app);
[mt, Nt] = size(test);
I = [];
testD = [];

if miss,
 [I J]=find(app==miss);
 [I2 J2]=find(test==miss);
end
completes=setdiff(1:ma,I);
donnees_continue=app(completes,continues);

% echantillonnage
[n,bornes,nbbornes,xx]=hist_ic(donnees_continue,critere);

% on re-distribue l'ensemble des donnees d'apprentissage continues
[n2,appD_continues]=histc_ic(app(:,continues),bornes);
if tt, [n2test,testD_continues]=histc_ic(test(:,continues),bornes); end

% on insere les donnees continues discretisees dans les matrices
appD=app;
if tt, testD=test; end
for i=1:length(continues)
    appD(:,continues(i))=appD_continues(:,i);
    if tt, testD(:,continues(i))=testD_continues(:,i); end
end

if miss,
 for k=1:length(I)
  app(I(k),J(k))=miss;
 end
 if tt, 
  for l=1:length(I2)
   testD(I2(l),J2(l))=miss;
  end
 end
end
appD = appD';
testD = testD';