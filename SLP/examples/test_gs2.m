% test Greedy Search with cache
clear all;

fprintf('\n=== Structure Learning with Greedy Search algorithm and cache implementation\n');

n=8;
names={ 'A' , 'S' , 'T' , 'L' , 'B' , 'O' , 'X' , 'D' };
carre=ones(1,n);
node_type={'tabular','tabular','tabular','tabular','tabular','tabular','tabular','tabular'};

m=6000;
bnet=mk_asia2_bnet;
data = cell(n,m);
for l = 1:m, data(:,l) = sample_bnet(bnet); end
data=cell2mat(data);
fprintf('Complete data have been created.\n\n');

seeddag = mk_rnd_dag(n);

fprintf('\t- Greedy Search (with cache)\n');
tmp=cputime;
L=300;
cache=score_init_cache(n,L);
[dag1, best_score, cache] =learn_struct_gs2(data,bnet.node_sizes,seeddag,'cache',cache,'scoring_fn','bayesian');
%cache
tmp=cputime-tmp;
fprintf('\t- Execution time : %3.2f seconds\n',tmp);

fprintf('\t- Greedy Search (without cache)\n');
tmp=cputime;
dag2=learn_struct_gs2(data,bnet.node_sizes,seeddag,'scoring_fn','bayesian');
tmp=cputime-tmp;
fprintf('\t- Execution time : %3.2f seconds\n',tmp);

fprintf('\t- Greedy Search (with cache and MWST initialisation)\n');
tmp=cputime;
Tdag=learn_struct_mwst(data, ones(n,1), bnet.node_sizes, node_type,'mutual_info',ceil(n*rand));
dag3=learn_struct_gs2(data,bnet.node_sizes,full(Tdag),'cache',cache,'scoring_fn','bayesian');
tmp=cputime-tmp;
fprintf('\t- Execution time : %3.2f seconds\n',tmp);


figure;[xx yy] = make_layout(bnet.dag);
yy=(yy-0.2)*.8/.6+.1;
xx=(xx-0.2833)*.8/.517+.1;
subplot(1,4,1), [xx yy]=draw_graph(bnet.dag,names,carre,xx,yy); %,carre);
title('ASIA original graph');
subplot(1,4,2), draw_graph(dag1,names,carre,xx,yy); %,carre);
title('GS (cache)');
subplot(1,4,3), [xx yy]=draw_graph(dag2,names,carre,xx,yy); %,carre);
title('GS (without cache)');
subplot(1,4,4), [xx yy]=draw_graph(dag3,names,carre,xx,yy); %,carre);
title('GS (with cache and MWST init)');
drawnow;
