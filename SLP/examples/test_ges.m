% test Greedy Equivalence Search with cache

clear all;

%add_BNT_to_path;

fprintf('\n=== Structure Learning with Greedy Equivalence Search algorithm and cache implementation\n');


bnet=mk_asia_bnet ;

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

seeddag = zeros(n,n);

fprintf('\t- Greedy Search (with cache)\n');
tmp=cputime;
L=500;
cache=score_init_cache(n,L);
[dag1, best_score, cache] =learn_struct_ges(data,bnet.node_sizes,seeddag,'cache',cache,'scoring_fn','bic');
%cache
tmp=cputime-tmp;
fprintf('\t- Execution time : %3.2f seconds\n',tmp);


figure;[xx yy] = make_layout(bnet.dag);
yy=(yy-0.2)*.8/.6+.1;
xx=(xx-0.2833)*.8/.517+.1;
subplot(1,3,1), [xx yy]=draw_graph(bnet.dag,names,carre,xx,yy); %,carre);
title('ASIA original graph');
subplot(1,3,2), draw_graph(dag1,names,carre,xx,yy); %,carre);
title('GES CPDAG (cache)');
subplot(1,3,3), draw_graph(cpdag_to_dag(dag1),names,carre,xx,yy); %,carre);
title('GES DAG (cache)');
drawnow;