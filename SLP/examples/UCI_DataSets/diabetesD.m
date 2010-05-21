load -ascii diabetes

[N, m] = size(diabetes)

class = N

%rand('state',0); randn('state',0);
%diabetes = diabetes(:,randperm(m));

Napp = 400
Ntest = m-Napp

app  = diabetes(:,1:Napp);size(app)
test = diabetes(:,Napp+1:end);size(test)

unique(app(class,:))
unique(test(class,:))

ns = max(diabetes')
clear diabetes

% N, ns(class), Napp, Ntest, mean(ns),
