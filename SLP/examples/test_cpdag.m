% Test Dag2Cpdag and Cpdag2Dag functions ...

close all;
clear all;

bnet=mk_asia_bnet ;

n=8;
names={ 'A' , 'S' , 'T' , 'L' , 'B' , 'O' , 'X' , 'D' };
carre=ones(1,n);

[xx yy] = make_layout(bnet.dag);
yy=(yy-0.2)*.8/.6+.1;
xx=(xx-0.2833)*.8/.517+.1;

cpdag=dag_to_cpdag(bnet.dag) ;
dag1 = cpdag_to_dag(cpdag) ;

figure; 
subplot(1,3,1), draw_graph(bnet.dag,names,carre,xx,yy); title('original ASIA');
subplot(1,3,2), draw_graph(cpdag,names,carre,xx,yy); title('cpasia=DAGtoCPDAG(asia)');
subplot(1,3,3), draw_graph(dag1,names,carre,xx,yy); title('CPDAGtoDAG(cpasia)');