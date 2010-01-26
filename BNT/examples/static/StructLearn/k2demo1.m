N = 4;
dag = zeros(N,N);
%C = 1; S = 2; R = 3; W = 4;
C = 4; S = 2; R = 3; W = 1; % arbitrary order
dag(C,[R S]) = 1;
dag(R,W) = 1;
dag(S,W)=1;

false = 1; true = 2;
ns = 2*ones(1,N); % binary nodes

bnet = mk_bnet(dag, ns);
bnet.CPD{C} = tabular_CPD(bnet, C, [0.5 0.5]);
bnet.CPD{R} = tabular_CPD(bnet, R, [0.8 0.2 0.2 0.8]);
bnet.CPD{S} = tabular_CPD(bnet, S, [0.5 0.9 0.5 0.1]);
bnet.CPD{W} = tabular_CPD(bnet, W, [1 0.1 0.1 0.01 0 0.9 0.9 0.99]);

seed = 0;
rand('state', seed);
randn('state', seed);
ncases = 100;
data = zeros(N, ncases);
for m=1:ncases
  data(:,m) = cell2num(sample_bnet(bnet));
end

order = [C S R W];
max_fan_in = 2;

%dag2 = learn_struct_K2(data, ns, order, 'max_fan_in', max_fan_in, 'verbose', 'yes');
  
sz = 5:5:50;
for i=1:length(sz)
  dag2 = learn_struct_K2(data(:,1:sz(i)), ns, order, 'max_fan_in', max_fan_in);
  correct(i) = isequal(dag, dag2);
end
correct

for i=1:length(sz)
  dag3 = learn_struct_K2(data(:,1:sz(i)), ns, order, 'max_fan_in', max_fan_in, 'scoring_fn', 'bic', 'params', []);
  correct(i) = isequal(dag, dag3);
end
correct


