load -ascii abalone

[N, m] = size(abalone)

class = N

%rand('state',0); randn('state',0);
%abalone = abalone(:,randperm(m));

Napp = 3133;
Ntest = m-Napp

app  = abalone(:,1:Napp);size(app)
test = abalone(:,Napp+1:end);size(test)

unique(app(class,:))
unique(test(class,:))

ns = max(abalone')
clear abalone

% N, ns(class), Napp, Ntest, mean(ns),
