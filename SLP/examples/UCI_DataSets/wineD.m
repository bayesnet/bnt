load -ascii wine

[N, m] = size(wine)

class = 1

rand('state',0); randn('state',0);
wine = wine(:,randperm(m));

Napp = 80
Ntest = m-Napp

app  = wine(:,1:Napp);size(app)
test = wine(:,Napp+1:151);size(test)

%unique(app(class,:))
%unique(test(class,:))

%app  = wine;
%test = wine;
%disp('test en apprentissage');

ns = max(wine')
clear wine

% N, ns(class), Napp, Ntest, mean(ns),
