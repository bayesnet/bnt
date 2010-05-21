load -ascii SPECT_A
load -ascii SPECT_T

[N, m] = size(SPECT_A);

class = N

app  = SPECT_A;size(app)
test = SPECT_T;size(test)

Napp = size(app,2);
Ntest = size(test,2);

unique(app(class,:))
unique(test(class,:))

ns1 = max(SPECT_A');
ns2 = max(SPECT_T');
ns = max(ns1, ns2)
clear SPECT_A SPECT_T ns1 ns2

% N, ns(class), Napp, Ntest, mean(ns),
