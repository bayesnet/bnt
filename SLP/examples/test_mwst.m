% test MWST

clear all;

fprintf('\n=== Structure Learning with Maximum Weight Spanning Tree\n');

bnet=mk_asia2_bnet ;

n=8;
names={ 'A' , 'S' , 'T' , 'L' , 'B' , 'O' , 'X' , 'D' };
carre=ones(1,n);
node_type={'tabular','tabular','tabular','tabular','tabular','tabular','tabular','tabular'};

m=10000;
bnet=mk_asia2_bnet;
data = cell(n,m);
for l = 1:m, data(:,l) = sample_bnet(bnet); end
data=cell2mat(data);
fprintf('Complete data have been created.');

tmp=cputime;
dag=learn_struct_mwst(data, ones(n,1), bnet.node_sizes, node_type,'mutual_info');
tmp=cputime-tmp;
fprintf('\t- Execution Time : %3.2f seconds\n\n',tmp);

figure; [xx yy] = make_layout(bnet.dag);
yy=(yy-0.2)*.8/.6+.1;
xx=(xx-0.2833)*.8/.517+.1;
subplot(1,2,1), [xx yy]=draw_graph(bnet.dag,names,carre,xx,yy); %,carre);
title('ASIA original graph');
subplot(1,2,2), draw_graph(dag,names,carre,xx,yy); %,carre);
title('MWST DAG');
drawnow;
