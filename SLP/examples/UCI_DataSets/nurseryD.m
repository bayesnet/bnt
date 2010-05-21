load -ascii nursery

[N, m] = size(nursery);

class = N

rand('state',0); randn('state',0);
nursery = nursery(:,randperm(m));

Napp = 8500
Ntest = m-Napp

app  = nursery(:,1:Napp);size(app)
test = nursery(:,Napp+1:end);size(test)

unique(app(class,:))
unique(test(class,:))

ns = max(nursery')
clear nursery

% N, ns(class), Napp, Ntest, mean(ns),
