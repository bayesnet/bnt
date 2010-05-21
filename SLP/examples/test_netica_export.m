%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Example 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=3;
dag=zeros(N);
dag([1 3],2)=1;

names = {'A1','B2','C3'};

bnet = mk_bnet(dag, [2 3 4], 'names', names);
%bnet = mk_bnet(dag, [2 3 4]);

bnet.CPD{1} = tabular_CPD (bnet, 1, [0.2 0.8] );
bnet.CPD{3} = tabular_CPD (bnet, 3, [0.1 0.2 0.3 0.4] );

      % 1 3 2 prob
tab = [[1 1 1 .2];...
       [2 1 1 .4];...
       [1 2 1 .1];...
       [2 2 1 .7];...
       [1 3 1 .5];...
       [2 3 1 .1];...
       [1 4 1 .75];...
       [2 4 1 .1];...
...
       [1 1 2 .3];...
       [2 1 2 .5];...
       [1 2 2 .1];...
       [2 2 2 .2];...
       [1 3 2 .25];...
       [2 3 2 .9];...
       [1 4 2 .10];...
       [2 4 2 .25];...
...
       [1 1 3 .5];...
       [2 1 3 .1];...
       [1 2 3 .8];...
       [2 2 3 .1];...
       [1 3 3 .25];...
       [2 3 3 .0];...
       [1 4 3 .15];...
       [2 4 3 .65]];

bnet.CPD{2} = tabular_CPD (bnet, 2, tab(:,4));

export_dnet(bnet,'test1', 1)

m = 20;
data = cell(N,m);
for l = 1:m, data(:,l) = sample_bnet(bnet); end
remov = rand(1,N*m)>.8;
for i=1:length(remov); if remov(i), data{i}=[];end; end

export_cases(data, names, 'test1', -9999)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Example 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=4;
dag=zeros(N);
dag([1 2 3],4)=1;
dag(1,2)=1;
ns = [2 2 2 2];

names = {'A1','B2','C3','D4'};

bnet = mk_bnet(dag, ns, 'names', names);

bnet.CPD{1} = tabular_CPD (bnet, 1, [0.2 0.8] );
bnet.CPD{2} = tabular_CPD (bnet, 2, [0.1 0.4 0.9 0.6] );
bnet.CPD{3} = tabular_CPD (bnet, 3, [0.5 0.5] );
bnet.CPD{4} = tabular_CPD (bnet, 4, [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2] );

export_dnet(bnet,'test2', 1)

m = 20;
data = cell(N,m);
for l = 1:m, data(:,l) = sample_bnet(bnet); end
remov = rand(1,N*m)>.8;
for i=1:length(remov); if remov(i), data{i}=[];end; end

export_cases(data, names, 'test2', -9999)

