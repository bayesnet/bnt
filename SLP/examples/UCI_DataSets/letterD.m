load -ascii letter_A
load -ascii letter_T

[N, m] = size(letter_A)

class = N

app  = letter_A;size(app)
test = letter_T;size(test)

Napp = size(app,2);
Ntest = size(test,2);

unique(app(class,:))
unique(test(class,:))

ns1 = max(letter_A');
ns2 = max(letter_T');
ns = max(ns1,ns2)
clear letter_A letter_T ns1 ns2

% N, ns(class), Napp, Ntest, mean(ns),
