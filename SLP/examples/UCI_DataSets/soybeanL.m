load -ascii soybeanD
load -ascii soybeanTD

% soybeanD(2:36,:) = soybeanD(2:36,:)+1;
% soybeanTD(2:36,:) = soybeanTD(2:36,:)+1;
% soybeanD(find(soybeanD==-9998))=-9999;
% soybeanTD(find(soybeanTD==-9998))=-9999;

[N, Napp] = size(soybeanD);
[N, Ntest] = size(soybeanTD);

N, m=Napp+Ntest, Napp, Ntest, 
class = 1

%rand('state',0); randn('state',0);
%abalone = abalone(:,randperm(m));

app  = soybeanD;
test = soybeanTD;

unique(app(class,:))
unique(test(class,:))

ns = max(soybeanD')
clear soybeanD soybeanTD

% N, ns(class), Napp, Ntest, mean(ns),
