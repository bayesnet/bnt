load -ascii house.dat
houseD=house';
clear house

[N, m] = size(houseD)

class = 1

%rand('state',0); randn('state',0);
%houseD = houseD(:,randperm(m));

Napp = ceil(m*2/3);
Ntest = m-Napp

app  = houseD(:,1:Napp);size(app)
test = houseD(:,Napp+1:end);size(test)

unique(app(class,:))
unique(test(class,:))

ns = max(houseD')
clear houseD
