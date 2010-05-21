load -ascii heart

[N, m] = size(heart)

class = N

%rand('state',0); randn('state',0);
%heart = heart(:,randperm(m));

Napp = 150
Ntest = m-Napp

app  = heart(:,1:Napp);size(app)
test = heart(:,Napp+1:end);size(test)

unique(app(class,:))
unique(test(class,:))

ns = max(heart')
clear heart

% N, ns(class), Napp, Ntest, mean(ns),
