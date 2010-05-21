function dag = learn_struct_tan_EM(data, class, node_sizes, root, prior, nbloopmax, thresh)
% LEARN_STRUCT_TAN_EM 
% dag = learn_struct_tan_EM(data, class, node_sizes, root, prior, nbloopmax, thresh)
%
% Learn TAN classifier for discrete variables from incomplete dataset
%
% Input : 
%   data{i,m} a cell where the node i in the case m,
%   class is the number of the class node,
%   node_sizes = 1 if gaussian node,                 (max on complete samples)
%   root is the futur root-node of the tree T.       (random)
%   prior = 1 to use uniform Dirichlet prior         (0)
%   nbloopmax = max loop number                      (ceil(log(N*log(N))))
%   thresh = the convergence test's threshold        (1e-3)
%
% Output :
%	bnet = the output bayesian network
%       Ebic = the espected BIC score of bnet given the data
%
% francois.olivier.c.h@gmail.com

[N, m]=size(data);
log_m = log(m);
if nargin<7, thresh = 1e-4; end
max_iter = 10; % for learn_struct_params
if nargin<6, nbloopmax = 15, end
if nargin<5, prior = 0; end
if nargin<4, root = ceil(N*rand(1)), end
if nargin<3,
  misv = -9999;
  data_mat = bnt_to_mat(data,misv);
  node_sizes = max(data_mat'),
  clear data_mat
end

discrete = 1:N;
nbloop = 0;
variab = mysetdiff(1:N, class);

[Bbest1, Sbest, Obest] = learn_struct_mwst_EM4tan(data(variab, :), 1:N-1, node_sizes(variab), prior, nbloopmax, thresh);
dag = zeros(N);
dag(variab, variab) = Bbest1.dag;
dag(class, variab) = 1;

%Bbest = mk_bnet(dag, discrete);
%for i=discrete, Bbest.CPD{i} = tabular_CPD(Bbest, i, 'prior_type', 'dirichlet', 'dirichlet_type', 'unif'); end
%engine = jtree_inf_engine(Bbest);
%Bbest = learn_params_em(engine, data, max_iter, thresh);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Bbest, Sbest, Obest] = learn_struct_mwst_EM4tan(data, discrete, node_sizes, prior, nbloopmax, thresh)

%fprintf('-- INITIALIZATION\n');
%%%%%%%%%%%%%
[N, m]=size(data);
 rand('state',sum(100*clock))

if nargin<6, thresh = 1e-4; end
if nargin<5, nbloopmax = 15; end
if nargin<4, prior = 0; end
if nargin<3, 
    misv = -9999;
    data_mat = bnt_to_mat(data,misv);
    node_sizes = max(data_mat'), 
end
if nargin<2, discrete = ones(1,N); end

nbloop = 1; last=0;
max_iter = 15; % for learn_struct_params

%%%%%%%%%%%%%
% Random Chain like DAG
T = diag(ones(N-1,1),1); T=T+T';
order=randperm(N); root=randperm(N); root=root(1);
[tmp order2]=sort(order); order2;
T2 = full(mk_rooted_tree(T,order2(root)));
torder=topological_sort(T2(order2,order2));
[tmp torder2]=sort(torder); torder2;

bnet0 = mk_bnet(T2(order2(torder),order2(torder)), node_sizes(torder));
for i=1:N, bnet0.CPD{i} = tabular_CPD(bnet0, i, 'prior_type', 'dirichlet', 'dirichlet_type', 'unif'); end % a priori -> change les espected counts

%%%%%%%%%%%%%
Sbest = -Inf; Bbest = bnet0; fini=0;

while not(fini)
    %fprintf('Loop %d : learning parameters...\n',nbloop);
    engine0=jtree_inf_engine(bnet0);  
    [bnet1, LL1, engine1] = learn_params_em(engine0, data(torder,:), max_iter, thresh);

    BIC0=0;
    for i=1:N, 
        xxx=struct(bnet1.CPD{i}); 
        BIC0=BIC0+bic_score_family(xxx.counts, xxx.CPT, xxx.nsamples);
    end

    if BIC0 < Sbest+ thresh*abs(Sbest),
        fini=1; 
    else
        Sbest = BIC0; 
        Bbest = bnet1;  
        Obest = torder2 ;
    %%%%%%%%%%%%%

    theta_Xi=cell(N,1);
    evidence = cell(1,N);
    [engine2, loglik] = enter_evidence(engine1, evidence);
    for j=1:N,
        SS= marginal_nodes (engine2,torder2(j)); 
        theta_Xi{j} = SS.T; 
    end

    theta_Xj_given_Xi = cell(N,N); 
    for i=1:N
        for vali = 1:node_sizes(i)
            evidence = cell(1,N);  evidence{torder2(i)} = vali;
            [engine2, loglik] = enter_evidence(engine1, evidence);
            for j=mysetdiff(1:N,i)
                SS= marginal_nodes (engine2,torder2(j)); 
                theta_Xj_given_Xi{j,i} = [theta_Xj_given_Xi{j,i}, SS.T];
            end
        end
    end

    BIC_mat=zeros(N,N);
    for i=1:N,
        BIC_mat(i,i)=bic_score_family(theta_Xi{i}*m,theta_Xi{i},m);
        for j=mysetdiff(1:N,i)
            theta_XjXi=(ones(node_sizes(i),1)*theta_Xi{j}').*theta_Xj_given_Xi{i,j};
            BIC_mat(i,j)= bic_score_family(m*theta_XjXi,theta_Xj_given_Xi{i,j},m);
        end
    end
    BIC_delta = BIC_mat-diag(BIC_mat)*ones(1,N);

    BIC1=0;
    for i = 1:N
        j = find(bnet0.dag(:,i)==1);
        if isempty(j)
            BIC1 = BIC1 + BIC_mat(torder2(i),torder2(i));
        else
            BIC1 = BIC1 + BIC_mat(torder2(i),torder2(j));
        end
    end

        %fprintf('     Creation of the new bnet                           ');tic
        %%%%%%%%%%%%%  
        T2 = minimum_spanning_tree(-BIC_delta);
        root=randperm(N); root=root(1);
        T3 = full(mk_rooted_tree(T2, root));
        torder=topological_sort(T3);
        [tmp torder2]=sort(torder); torder2;

        bnet0 = mk_bnet(T3(torder,torder), node_sizes(torder));
        %BIC=0;
        for i = 1:N
            j = find(T3(:,i)==1);
            if isempty(j)
                bnet0.CPD{torder2(i)} = tabular_CPD(bnet0, torder2(i), 'CPT', theta_Xi{i}, 'prior_type', 'dirichlet', 'dirichlet_type', 'unif');
                %BIC = BIC + BIC_mat(i,i);
            else
                bnet0.CPD{torder2(i)} = tabular_CPD(bnet0, torder2(i), 'CPT', theta_Xj_given_Xi{i,j}, 'prior_type', 'dirichlet', 'dirichlet_type', 'unif');
                %BIC=BIC + BIC_mat(i,j);
            end
        end
        nbloop=nbloop+1;
        fprintf('================================================================================\n');
        if nbloop>nbloopmax
            fini=1;
        end
    end
end
