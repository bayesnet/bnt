load -ascii tae

[N, m] = size(tae)

class = 6

%rand('state',0); randn('state',0);
%tae = tae(:,randperm(m));

Napp = 100
Ntest = m-Napp

app  = tae(:,1:Napp);size(app)
test = tae(:,Napp+1:151);size(test)

%unique(app(class,:))
%unique(test(class,:))

%app=tae;
%test=tae;
%disp('test en apprentissage');

ns = max(tae')
clear tae

% N, ns(class), Napp, Ntest, mean(ns),
