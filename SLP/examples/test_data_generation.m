% gener_mcardat

%%Jouet 1
N=5;
dag=zeros(N);
dag(1,2)=1;
dag(2,[3 4])=1;
dag(4,5)=1;
node_sizes= [2 3 4 5 2]; ns = node_sizes;
bnet = mk_bnet(dag, node_sizes);
bnet.CPD{1} = tabular_CPD(bnet, 1, [0.2 0.8]);
bnet.CPD{2} = tabular_CPD(bnet, 2, [0.1 0.6 0.3 0.9 0.4 0.7]);
bnet.CPD{3} = tabular_CPD(bnet, 3, [0.3 0.5 0.7 0.2   0.2 0.1 0.05 0.5   0.5 0.4 0.25 0.3]);
bnet.CPD{4} = tabular_CPD(bnet, 4, [0.35 0.05 0.65 0.15 0.25  0.15 0.15 0.10 0.25 0.45  0.5 0.8 0.25 0.6 0.3]);
bnet.CPD{5} = tabular_CPD(bnet, 5, [0.25 0.15  0.15 0.25  0.10 0.15  0.2 0.05  0.3 0.4]);
bnet_orig=bnet

base_proba=0.2;

%%%%%%%% MCAR

bnet_miss = gener_MCAR_net(bnet_orig, base_proba);

carre=zeros(1,3*N);
names={'X1','X2','X3','X4','X5','R1','R2','R3','R4','R5','M1','M2','M3','M4','M5'};
xx=[.1,.3,.5,.7,.9,.05,.25,.45,.65,.85,.15,.35,.55,.75,.95];
yy=[.85,.70,.90,.75,.80,.65,.40,.60,.45,.50,.2,.2,.2,.2,.2];
figure;draw_graph(bnet_miss.dag,names,carre,xx,yy)

[data, comp_data, bnet_miss, taux, bnet_orig, notok] = gener_data_from_bnet_miss(bnet_miss, 500, base_proba ,0,0);

OK = ~notok

%%%%%%%% MAR

bnet_miss = gener_MAR_net(bnet_orig, base_proba);

carre=zeros(1,3*N);
names={'X1','X2','X3','X4','X5','R1','R2','R3','R4','R5','M1','M2','M3','M4','M5'};
xx=[.1,.3,.5,.7,.9,.05,.25,.45,.65,.85,.15,.35,.55,.75,.95];
yy=[.85,.70,.90,.75,.80,.65,.40,.60,.45,.50,.2,.2,.2,.2,.2];
figure;draw_graph(bnet_miss.dag,names,carre,xx,yy)

[data, comp_data, bnet_miss, taux, bnet_orig, notok] = gener_data_from_bnet_miss(bnet_miss, 500, base_proba ,0,0);

OK = ~notok
