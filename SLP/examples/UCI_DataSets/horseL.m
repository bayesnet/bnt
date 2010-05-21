load -ascii horseD
load -ascii horseTD

% horse = horse';
% horse(2,:)=(horse(2,:)-1)/8+1;
% horse(27,:)=horse(27,:)/2209+1;
% horseT = horseT';
% horseT(2,:)=(horseT(2,:)-1)/8+1;
% horseT(27,:)=horseT(27,:)/2209+1;
% [horseD, horseTD] = discretization(3, [3,4,5,6,16,19,20,22,25], -9999, horse, horseT);
% horseD(26,:)= nominalval(horseD(26,:));
% horseTD(26,:)= nominalval(horseTD(26,:));

[N, Napp] = size(horseD);
ns=max([horseD,horseTD]')
UNI=find(ns==1);
horseD=horseD(mysetdiff(1:N,UNI),:);
horseTD=horseTD(mysetdiff(1:N,UNI),:);
ns=max([horseD,horseTD]')

[N, Napp] = size(horseD);
[N, Ntest] = size(horseTD);

N, m=Napp+Ntest, Napp, Ntest, 
class = 1

%rand('state',0); randn('state',0);

app  = horseD;
test = horseTD;
fotest = test(:,9); % la classe est absente :-/
newtest = app(:,299);
app(:,299) = fotest;
test(:,9) = newtest;

unique(app(class,:))
unique(test(class,:))

clear horseD horseTD fotest newtest

% N, ns(class), Napp, Ntest, mean(ns),
