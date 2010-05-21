load -ascii pen_A
load -ascii pen_T

[N, m] = size(pen_A);

class = N

app  = pen_A;size(app)
test = pen_T;size(test)

Napp = size(app,2);
Ntest = size(test,2);

unique(app(class,:))
unique(test(class,:))

ns1 = max(pen_A');
ns2 = max(pen_T');
ns = max(ns1, ns2)
clear pen_A pen_T ns1 ns2

% N, ns(class), Napp, Ntest, mean(ns),
