load -ascii mushrooms.dat
mushroomsD=mushrooms';
clear mushrooms

[N, m] = size(mushroomsD);
class = 1

%rand('state',0); randn('state',0);
%abalone = abalone(:,randperm(m));

for node = 1:N
  UNI=setdiff(unique(mushroomsD(node,:)),-9999);
  for val = 1:length(UNI)
    [I,J]=find(mushroomsD(node,:)==UNI(val));
    mushroomsD(node,J)=val;
  end
end

ns = max(mushroomsD');
seul=find(ns==1);
mushroomsD=mushroomsD(setdiff(1:N,seul),:);
[N, m] = size(mushroomsD)
ns = max(mushroomsD')

Napp = ceil(m*2/3);
Ntest = m-Napp

app  = mushroomsD(:,1:Napp);size(app)
test = mushroomsD(:,Napp+1:end);size(test)

unique(app(class,:))
unique(test(class,:))

clear mushroomsD seul UNI node I J
