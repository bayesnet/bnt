load -ascii german

[N, m] = size(german)

class = N

%rand('state',0); randn('state',0);
%german = german(:,randperm(m));

Napp = 600
Ntest = m-Napp

app  = german(:,1:Napp);size(app)
test = german(:,Napp+1:end);size(test)

unique(app(class,:))
unique(test(class,:))

ns = max(german')
clear german

% N, ns(class), Napp, Ntest, mean(ns),
