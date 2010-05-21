load -ascii segment

[N, m] = size(segment)

class = N

%rand('state',0); randn('state',0);
%segment = segment(:,randperm(m));

Napp = 1400
Ntest = m-Napp

app  = segment(:,1:Napp);size(app)
test = segment(:,Napp+1:end);size(test)

unique(app(class,:))
unique(test(class,:))

ns = max(segment')
clear segment

% N, ns(class), Napp, Ntest, mean(ns),
