% francois.olivier.c.h@gmail.com

%ddd = datestr(now);
%ddd([12 15 18])='-' ;
%fnd=[ddd '.txt'];
%diary(fnd)

%dbstop if error

clear all;
close all;

 rand('state',sum(100*clock))

 nbloopmax = 5;        % number of loop max in MWST-EM


names={ 'A' , 'S' , 'T' , 'L' , 'B' , 'O' , 'X' , 'D' };
node = struct('visit', 1, ...
    'smoking', 2, ...
    'tuberculosis', 3, ...
    'bronchitis', 5, ...
    'lung', 4, ...
    'ou', 6, ...
    'Xray', 7, ...
    'dyspnoea', 8);

adjacency = zeros(8);
adjacency([node.visit], node.tuberculosis) = 1;
adjacency([node.smoking], node.lung) = 1;
adjacency([node.lung node.tuberculosis], node.ou) = 1;
adjacency([node.ou], node.Xray) = 1;
adjacency([node.smoking], node.bronchitis) = 1;
adjacency([node.bronchitis node.ou], node.dyspnoea) = 1;
carre=ones(1,8);

  figure(1); [xx yy] = make_layout(adjacency);
    yy=(yy-0.2)*.8/.6+.1;
    xx=(xx-0.2833)*.8/.517+.1;
    subplot(2,2,1), [xx yy]=draw_graph(adjacency,names,carre,xx,yy); %,carre);
    title('ASIA net.');

fprintf('\n============================= Test MWST-EM\n');

n=8;
m=500;
bnet=mk_asia2_bnet;
data = cell(n,m);
for l = 1:m, data(:,l) = sample_bnet(bnet); end
asiab=cell2mat(data);
fprintf('Complete data have been created.');

 DM = 0.1;
 BD0 = asiab;
 node_sizes = max(BD0');
 [N, m]=size(BD0);
 rand('state',0); randn('state',0);
 vide = rand(size(BD0))<(1-DM);
 data=BD0.*vide;
 data = mat_to_bnt(data,0);

% N=4;
% dagO = diag(ones(N-1,1),1); dag0(1,3)=1;
% figure(1), subplot(4,4,1), title('theoritical'), draw_graph(dagO);
%
% node_sizes=2*ones(1,N);
discrete = ones(1,N);
%
% bnetO = mk_bnet(dagO, node_sizes);
% bnetO.CPD{1} = tabular_CPD(bnetO, 1, 'CPT',[0.2 0.8]);
% bnetO.CPD{2} = tabular_CPD(bnetO, 2, 'CPT',[0.4 0.7 0.6 0.3]);
% bnetO.CPD{3} = tabular_CPD(bnetO, 3);
% bnetO.CPD{4} = tabular_CPD(bnetO, 4, 'CPT', [0.5 0.8 0.5 0.2]);
%
% m = 1000; DM = 0.1;
%
% for l=1:m, dataO(:,l) = sample_bnet(bnetO); end
% rand('state',0); randn('state',0);
% vide = rand(size(dataO))<(1-DM);
% data = bnt_to_mat(dataO);
% data = data.*vide;
% data = mat_to_bnt(data,0);
% clear dataO;

fprintf('Missing data percentage : %3.1f%%\n',100*DM);

%     engine0=jtree_sparse_inf_engine(bnetO);
%     [bnet1, LL1, engine1] = learn_params_em(engine0, data);
% BIC0=0;
%     for i=1:N,
%         xxx=struct(bnet1.CPD{i});
%         BIC0=BIC0+bic_score_family(xxx.counts, xxx.CPT, xxx.nsamples);
%     end
%     fprintf('%5.2f\n',BIC0);

%root = 1;
prior = 0;
  tmp=cputime;

[BT_J11, Sbest0] = learn_struct_mwst_EM(data, discrete, node_sizes, prior, nbloopmax);
    tmp=cputime-tmp;
    fprintf('\tMWST-EM algorithm spent %3.2f secondes\n',tmp);

figure(1),  subplot(2,2,2), draw_graph(BT_J11.dag,names,carre,xx,yy); %,carre);
    title('MSWT-EM');

    fprintf('\n============================= Test AM-SEM\n');

    G0 = zeros(N,N);
    B0 = mk_bnet(G0, node_sizes);
    for i=1:N
        B0.CPD{i} = tabular_CPD(B0, i, 'prior_type', 'dirichlet', 'dirichlet_weight', 0);%1, 'dirichlet_type','BDeu');
    end

    tmp=cputime;
    max_loop = 10;

    [B0, order, best_score] = learn_struct_EM(B0, data, max_loop);
    G1 = B0.dag;
    [xxx oo]=sort(order);
    dag=G1(oo,oo);

    tmp=cputime-tmp;
    fprintf('\tSEM algorithm spent %3.2f secondes\n',tmp);


    subplot(2,2,3), draw_graph(dag,names,carre,xx,yy); %,carre);
    title('AMS-EM');

    fprintf('\n============================= Test AM-SEM+T\n');

    G0 = BT_J11.dag;
    B0 = mk_bnet(G0, node_sizes);
    for i=1:N
        B0.CPD{i} = tabular_CPD(B0, i, 'prior_type', 'dirichlet', 'dirichlet_weight', 0);%1, 'dirichlet_type','BDeu');
    end

    tmp=cputime;
    max_loop = 10;

    [B0, order, best_score] = learn_struct_EM(B0, data, max_loop);
    G1 = B0.dag;
    [xxx oo]=sort(order);
    dag=G1(oo,oo);

    tmp=cputime-tmp;
    fprintf('\tSEM+T algorithm spent %3.2f secondes\n',tmp);


    subplot(2,2,4), draw_graph(dag,names,carre,xx,yy); %,carre);
    title('AMS-EM+T');


%diary off
