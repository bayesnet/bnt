load -ascii hepatitis.dat
hepatitisD=hepatitis';
clear hepatitis

[N, m] = size(hepatitisD)
misv=-9999;

class = 1

%rand('state',0); randn('state',0);
%abalone = abalone(:,randperm(m));

ns = max(hepatitisD');
Cont=find(ns>2);
for node = Cont
  [I,J]=find(hepatitisD(node,:)~=misv);
  [N,EDGES,NBEDGES,XECHAN] = hist_ic(hepatitisD(node,J),3);
  hepatitisD(node,J)=XECHAN;
end

[N, m] = size(hepatitisD);
ns = max(hepatitisD')

Napp = 90;
Ntest = m-Napp

app  = hepatitisD(:,1:Napp);size(app)
test = hepatitisD(:,Napp+1:end);size(test)

unique(app(class,:))
unique(test(class,:))

clear hepatitisD node Cont EDGES NBEDGES XECHAN

% N, ns(class), Napp, Ntest, mean(ns),
