
clear all
B0 = 1; Rtriple = 2; Damnio = 3;
B1 = 4; Ramnio = 5; Dabort = 6; 
B2 = 7; U = 8;

N = 8;
dag = zeros(N,N);
dag(B0, [Rtriple B1 Ramnio]) = 1;
dag(Rtriple, [Damnio Dabort]) = 1;
dag(Damnio, [B1 Ramnio]) = 1;
dag(B1, B2) = 1;
dag(Ramnio, [Dabort U]) = 1;
dag(Dabort, B2) = 1;
dag(B2, U) = 1;



ns = zeros(1,N);
ns(B0) = 2;
ns(B1) = 3;
ns(B2) = 4;
ns(Rtriple) = 2;
ns(Ramnio) = 3;
ns(Damnio) = 2;
ns(Dabort) = 2;
ns(U) = 1;

limid = mk_limid(dag, ns, 'chance', [B0 B1 B2], ...
		 'decision', [Damnio Dabort], 'utility', [U]);

% states of nature
healthy = 1; downs = 2; miscarry = 3; aborted = 4;
% test results
pos = 1; neg = 2; unk = 3;
% actions
yes = 1; no = 2;

% Prior probability baby has downs syndrome
tbl = zeros(2,1);
p = 1/1000; % from www.downs-syndrome.org.uk figure
p = 24/10000; % www-personal.umich.edu/~bobwolfe/560/review/Downs.pdf (for women agen 35-40)
tbl(healthy) = 1-p;
tbl(downs) = p;
limid.CPD{B0} = tabular_CPD(limid, B0, tbl);

% Reliability of triple screen test
% Unreliable sensor
% B0 -> Rtriple
tbl = zeros(2,2); % Rtriple = pos, neg
p = 0.5;  % high false positive rate (guess)
tbl(healthy, :) = [p 1-p];
p = 0.6; % low detection rate (march of dimes figure)
tbl(downs, :) = [p 1-p]; 
limid.CPD{Rtriple} = tabular_CPD(limid, Rtriple, tbl);

limid.CPD{Damnio} = tabular_decision_node(limid, Damnio);

% Effect of amnio on baby  B0,Damnio -> B1
 % 1/200 risk of miscarry 
p = 1/200; % (march of dimes figure)
tbl = zeros(2, 2, 3); % B1 = healthy, downs, miscarry
tbl(healthy, no, :) =  [1     0     0];
tbl(downs, no, :) =    [0     1     0];
tbl(healthy, yes, :) = [1-p     0   p];
tbl(downs, yes, :) =   [0     1-p   p];
limid.CPD{B1} = tabular_CPD(limid, B1, tbl);

% Reliability of amnio  B0, Damnio -> Ramnio
% Perfect sensor
tbl = zeros(2,2,3); % Ramnio = pos, neg, unk
tbl(:, no, :) =        repmat([0 0 1], 2 ,1);
tbl(healthy, yes, :) = [0 1 0]; 
tbl(downs, yes, :) =   [1 0 0]; 
limid.CPD{Ramnio} = tabular_CPD(limid, Ramnio, tbl);

limid.CPD{Dabort} = tabular_decision_node(limid, Dabort);

% Effect of abortion on baby  B1, Dabort -> B2
tbl = zeros(3, 2, 4); % B2 = healthy, downs, miscarry, aborted
tbl(:, yes, :) =       repmat([0 0 0 1], 3, 1);
tbl(healthy, no, :) =  [1 0 0 0];
tbl(downs, no, :) =    [0 1 0 0];
tbl(miscarry, no, :) = [0 0 1 0];
limid.CPD{B2} = tabular_CPD(limid, B2, tbl);

% Utility U(Ramnio, B2)
tbl = zeros(3, 4);
tbl(:, healthy) = 5000;
tbl(:, downs) = -50000;
tbl(:, miscarry) = -1000;
tbl(:, aborted) = -1000;

if 0
%tbl(unk, miscarry) = 0; % this case is impossible
tbl(pos, miscarry) = -1;
tbl(neg, miscarry) = -1000;
if 1
  tbl(unk, aborted) = -100;
  tbl(pos, aborted) = -1;
  tbl(neg, aborted) = -500;
else % pro-life utility fn
  tbl(unk, aborted) = -500000;
  tbl(pos, aborted) = -500000;
  tbl(neg, aborted) = -500000;
end 
end

limid.CPD{U} = tabular_utility_node(limid, U,  tbl);



engine = jtree_limid_inf_engine(limid);
[strategy, MEU] = solve_limid(engine);

% Rtriple U(Damnio=1=yes)  U(Damnio=2=no)
% 1=pos    0               1
% 2=neg    0               1
dispcpt(strategy{Damnio})
if isequal(strategy{Damnio}(1,:), strategy{Damnio}(2,:))
  % Rtriple result irrelevant
  doAmnio = argmax(strategy{Damnio}(1,:))
else
  doAmnio = 1;
end

% Rtriple Ramnio U(Dabort=yes=1) U(Dabort=no=2)
% 1=pos   1=pos  1               0
% 2=neg   1=pos  1               0
% 1=pos   2=neg  0               1
% 2=neg   2=neg  0               1
% 1=pos   3=unk  0               1
% 2=neg   3=unk  0               1
dispcpt(strategy{Dabort})

