load -ascii contrasep

[N, m] = size(contrasep);

class = N

%rand('state',0); randn('state',0);
%contrasep = contrasep(:,randperm(m));

Napp = 850
Ntest = m-Napp

app  = contrasep(:,1:Napp);size(app)
test = contrasep(:,Napp+1:end);size(test)

unique(app(class,:))
unique(test(class,:))

ns = max(contrasep')
clear contrasep

% N, ns(class), Napp, Ntest, mean(ns),
