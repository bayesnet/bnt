load -ascii australian

[N, m] = size(australian)

class = 15

%rand('state',0); randn('state',0);
%australian = australian(:,randperm(m));

Napp = 400
Ntest = m-Napp

app  = australian(:,1:Napp);size(app)
test = australian(:,Napp+1:end);size(test)

unique(app(class,:))
unique(test(class,:))

ns = max(australian')
clear australian

% N, ns(class), Napp, Ntest, mean(ns),
