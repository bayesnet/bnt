load -ascii car

[N, m] = size(car);

class = N

rand('state',0); randn('state',0);
car = car(:,randperm(m));

Napp = 1152
Ntest = m-Napp

app  = car(:,1:Napp);size(app)
test = car(:,Napp+1:end);size(test)

unique(app(class,:))
unique(test(class,:))

ns = max(car')
clear car

% N, ns(class), Napp, Ntest, mean(ns),
