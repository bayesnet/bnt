load -ascii zoo

[N, m] = size(zoo)

class = 17

rand('state',0); randn('state',0);
zoo = zoo(:,randperm(m));

Napp = 60
Ntest = m-Napp

app  = zoo(:,1:Napp);size(app)
test = zoo(:,Napp+1:101);size(test)

%unique(app(class,:))
%unique(test(class,:))

%app  = zoo;
%test = zoo;
%disp('Test en apprentissage');

ns = max(zoo')
clear zoo

% N, ns(class), Napp, Ntest, mean(ns),
