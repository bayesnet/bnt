load thyroid_app
load thyroid_test
thyroid_test=thyroid_test';

[N, Napp] = size(thyroid_app);
[N, Ntest] = size(thyroid_test);

N, m=Napp+Ntest, Napp, Ntest, 
class = 1

%rand('state',0); randn('state',0);
%abalone = abalone(:,randperm(m));

app  = thyroid_app;
test = thyroid_test;

unique(app(class,:))
unique(test(class,:))

ns = max(thyroid_app')
clear thyroid_app thyroid_test
