load -ascii monks_A2
load -ascii monks_T

[N, m] = size(monks_T);

class = 1

app  = monks_A2;size(app)

test = monks_T;size(test)

Napp = size(app,2);
Ntest = size(test,2);

unique(app(class,:))
unique(test(class,:))

ns = max(monks_T');
clear monks_A2 monks_T

% N, ns(class), Napp, Ntest, mean(ns),
