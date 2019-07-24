% test PC + Chi2 statistical test

clear all;

test='pearson'; % Pearson's Chi2
%test='LRT';     % Likelihood ratio test G2
fprintf('\n=== Structure Learning with PC algorithm (and %s CHI2 test)\n',test);

bnet=mk_asia2_bnet ;

n=8;
names={ 'A' , 'S' , 'T' , 'L' , 'B' , 'O' , 'X' , 'D' };
carre=ones(1,n);
node_type={'tabular','tabular','tabular','tabular','tabular','tabular','tabular','tabular'};

m=2000;
bnet=mk_asia2_bnet;
data = cell(n,m);
for l = 1:m, data(:,l) = sample_bnet(bnet); end
data=cell2mat(data);
fprintf('Complete data have been created.');

tmp=cputime;
PC.dag = learn_struct_pdag_pc('cond_indep_chisquare',n,n-2,data,test);
tmp=cputime-tmp;
fprintf('\t- Execution time : %3.2f seconds\n',tmp);

figure;
[xx yy] = make_layout(bnet.dag);
yy=(yy-0.2)*.8/.6+.1;
xx=(xx-0.2833)*.8/.517+.1;
subplot(1,3,1), [xx yy]=draw_graph(bnet.dag,names,carre,xx,yy); %,carre);
title('ASIA original graph');
subplot(1,3,2), draw_graph(abs(PC.dag),names,carre,xx,yy); %,carre);
title('PC PDAG');
dag = cpdag_to_dag(abs(PC.dag));
subplot(1,3,3), draw_graph(dag,names,carre,xx,yy);
title('PC DAG');
drawnow;