function [dag,dag_score,confusion_matrix,correct_rate] ...
    = learn_struct_and_classified(data,classify_node_num, struct_algorithm,...
    classify_algorithm,varargin)
% Input:
%     data : data(r,c) = value of node r in case c (can be a cell array);
%
%     classify_node_num : classify node's num;
%
%     struct_algorithm  : which struct algorithm you use to learn struct,we
%                         support total 15 struct algorithm,like PC, IC,
%                         BNPC, MWST, NAVIE, TAN, FAN, K2, MCMC, GS2, GES,
%                         MWST-K2, HC, GS+T, MWST_EM, TAN_EM struct
%                         algorithm;
%
%     classify_algorithm : which classify algorithm you use to classified,
%                          we support 3 classify algorithm,like HOLD_OUT,
%                          CV-5, CV-10;
%
%     varargin :
%         node_flag  :  default 'A';
%         label      : default {},every node significance;
%         score      : in some struct algorithm, we need score ,like 'bic',
%                     'aic','mutual_info','bayesian' and so on. And when we
%                      get a dag's socre, we need set a score method.default
%                     'bic';
%         root       : in some struct algorithm, we need set a root ,like
%                     'TAN',default root is equal 1;
%       max_fan_in   : in some struct algorithm, we need set max_fan_in, like
%                    'K2',it means the largest number of parents we allow
%                     per node [N];
%       order        : in 'K2' algorithm, we need a order. order(i) is the
%                     i'th node in the topological ordering, if not
%                     specified order. we use random order in 'K2' algorithm;
%       scoring_mwst : in 'MWST-K2' algorithm, score in mwst,default 'bic';
%       scoring_K2   : in 'MWST-K2' algorithm, score in K2, default 'bic';
%       mwst_k2_order: in 'MWST-K2' algorithm, order has two type "+" or
%                      '-',default "+";
%       nsamples     : in 'MCMC' algorithm, we set nsamples, it means number
%                      of samples to draw from the chain after burn-in.
%                      default round(length(data)/5);
%       prior        : in some 'XX_EM' algorithm, we set prior.when it
%                       equals 1 to use uniform Dirichlet prior ,default 0;
%       nbloopmax    : in some 'XX_EM' algorithm, max loop number,
%                       default ceil(log(N*log(N)));
%       thresh       : in some 'XX_EM' algorithm, the convergence test's
%                      threshold, default 1e-3;
%       classify_times : classify times,default 10;
% Output:
%     dag               :  the dag by learn struct;
%     dag_score         : dag 's score,default 'bic' score;
%     confusion_matrix  : confusion matrix
%     correct_rate      : correct rate
%
% For example
%    pass in this moment
%
% other
%      I love math,I love Candy too.My dog who named Candy went
%      to heaven when I programming, I miss Candy so much.
% Written By WANGXin(growlithe1205@gmail.com)
%

% set default parameters
node_flag = 'A';
label = {};
score = 'bic';
% get variable parameters
if ~isempty(varargin)
    args = varargin;
    nargs = length(args);
    
    if length(args) > 0
        if ischar(args{1})
            for i=1:2:nargs
                switch args{i}
                    case 'node_flag'
                        node_flag = args{i+1};
                    case 'label'
                        label = args{i+1};
                    case 'score'
                        score = args{i+1};
                    case 'params'
                        if isempty(args{i+1})
                            params = cell(1,n);
                        else
                            params = args{i+1};
                        end
                end
            end
        end
    end
end

% data_check
[r,~] = check_data_integer(data);

% date_preprocess
[~, node_sizes, node_type] = data_process(r,data,node_flag);

% set start time
start_time = cputime;

% learn_struct
dag = learn_struct(data,struct_algorithm, node_sizes, node_type,...
    classify_node_num,varargin);

% get learn struct time
learn_struct_time = cputime - start_time;

% classify
[confusion_matrix,correct_rate] = classified_test(dag,classify_node_num,data,...
    classify_algorithm,varargin);

% get classifiy time
classified_time = cputime - learn_struct_time - start_time;

% get dag score
dag_score = score_dags(data, node_sizes, {dag},'scoring_fn',score);

% show
fprintf('learn struct time is %10.5f s \n',learn_struct_time);
fprintf('classified  time  is %10.5f s \n',classified_time);
fprintf('dag score in  %s is %10.5e \n', score, dag_score);

if isempty(label)
    draw_graph(dag);
else
    [~] = node_dependent(dag,label);
    draw_graph(dag,label);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check data
function [r,c] = check_data_integer(data)

[r,c] = size(data);
if (r == 1)
    disp('The row of data can not equal 1 \n');
    return
end

for i=1:r
    for j=1:c
        if(data(i,j) ~= fix(data(i,j)))
            fprintf('[%d,%d] is not integer \n',i,j)
            break
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data process
function [bnet ,node_sizes, node_type] = data_process(r,...
    data, node_flag)

dag = zeros(r,r);
[node_sizes,node_name,node_type] = get_node(data,node_flag);
bnet = mk_bnet(dag, node_sizes, 'names', node_name, 'discrete', 1:r);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get node sizes
function [node_sizes,node_name,node_type] = get_node(data,node_flag)

[r,~] = size(data);

node_class_num = cell(1,r);
node_sizes = cell(1,r);
node_name = cell(1,r);
node_type = cell(1,r);
for i = 1:r
    eval([node_flag, num2str(i) '=num2str(i);']);
    node_name{i} = [node_flag, num2str(i)];
    node_class_num{i} = unique(data(i,:));
    node_sizes{i} = length(node_class_num{i});
    node_type{i} = 'tabular';
end

node_sizes = cell2mat(node_sizes);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% learn struct
function dag = learn_struct(data, struct_algorithm, node_sizes,...
    node_type,classify_node_num, varargin)

[r,~] = size(data);

% default parameters
score = 'bic';
root = 1;
if(root == classify_node_num)
    root = root + 1;
end
max_fan_in = r;
order = randperm(length(node_sizes));
nsamples = round(length(data)/5);
scoring_mwst = 'bic';
scoring_K2 = 'bic';
mwst_k2_order = '+';
prior = 0;
nbloopmax = 30;
thresh = 1e-4;
tr_method = 'auto';

% get parameters
if ~isempty(varargin{1,1})
    args = varargin{1,1}(1,:);
    nargs = length(args);
    if length(args) > 0
        if ischar(args{1})
            for i=1:2:nargs
                switch args{i}
                    case 'score'
                        score = args{i+1};
                    case 'root'
                        root = args{i+1};
                    case 'max_fan_in'
                        max_fan_in = args{i+1};
                    case 'order'
                        order = args{i+1};
                    case 'nsamples'
                        nsamples = args{i+1};
                    case 'scoring_mwst'
                        scoring_mwst = args{i+1};
                    case 'scoring_K2'
                        scoring_K2 = args{i+1};
                    case 'mwst_k2_order'
                        mwst_k2_order = args{i+1};
                    case 'prior'
                        prior = args{i+1};
                    case 'nbloopmax'
                        nbloopmax = args{i+1};
                    case 'thresh'
                        thresh = args{i+1};
                    case 'tr_method'
                        tr_method = args{i+1};
                    case 'params'
                        if isempty(args{i+1})
                            params = cell(1,n);
                        else
                            params = args{i+1};
                        end
                end
            end
        end
    end
end

switch struct_algorithm
    
    case 'PC'
        % pc algorithm
        PC.dag = learn_struct_pdag_pc('cond_indep_chisquare',r,r-2,data);
        dag = cpdag_to_dag(abs(PC.dag));
        
    case 'IC'
        % ic algorithm
        PC.dag = learn_struct_pdag_ic_star('cond_indep_chisquare',r,r-2,data);
        dag = cpdag_to_dag(abs(PC.dag));
        
    case 'BNPC'
        % bnpc algorithm
        dag = learn_struct_bnpc(data);
        
    case 'MWST'
        % mwst algorithm
        dag = learn_struct_mwst(data, ones(r,1), node_sizes, node_type,...
            score, classify_node_num);
        
    case 'NAIVE'
        % naive bayes algorithm
        dag = mk_naive_struct(r,classify_node_num);
        
    case 'TAN'
        % tan algorithm
        dag = learn_struct_tan(data, classify_node_num, root, node_sizes,...
            score);
        
    case 'FAN'
        % fan algorithm
        dag = learn_struct_fan(data, classify_node_num, node_sizes,...
            score,tr_method);
        
    case 'K2'
        % k2 algorithm order random
        dag = learn_struct_K2(data,node_sizes,order,max_fan_in,...
            'scoring_fn',score);
        
    case 'MCMC'
        % mcmc algorithm and get optimal dag
        dag = get_optimal_dag_by_mcmc(data,node_sizes,'nsamples',nsamples,...
            'scoring_fn',score);
        
    case 'GS2'
        % gs2 algorithm
        seeddag = full(learn_struct_mwst(data, ones(r,1), node_sizes, node_type));
        dag = learn_struct_gs2(data, node_sizes, seeddag);
        
    case 'GES'
        % ges algorithm
        dag = learn_struct_ges(data, node_sizes, score);
        
    case 'MWST-K2'
        % mwst_k2 algorithm
        dag = learn_struct_mwst_K2(data, node_sizes, node_type,classify_node_num,...
            'order',mwst_k2_order,'scoring_mwst',scoring_mwst,'scoring_K2',...
            scoring_K2,'max_fan_in',max_fan_in);
        
    case 'HC'
        % hc algorithm
        seeddag = full(learn_struct_mwst(data, ones(r,1), node_sizes, node_type));
        dag = learn_struct_hc(data,node_sizes,seeddag);
        
    case 'GS+T'
        % gs_t algorithm
        seeddag = full(learn_struct_mwst(data, ones(r,1), node_sizes, node_type));
        dag = learn_struct_gs2(data,node_sizes,seeddag);
        
    case 'MWST_EM'
        % mwst_em algorithm
        data = num2cell(data);
        bnet = learn_struct_mwst_EM(data,ones(r,1), node_sizes, prior,...
            nbloopmax,thresh);
        dag = bnet.dag;
        
    case 'TAN_EM'
        % tan_em algorithm
        data = num2cell(data);
        bnet = learn_struct_tan_EM(data, classify_node_num, node_sizes, root,...
            prior, nbloopmax, thresh);
        dag = bnet;
        
    otherwise
        fprintf('%s struct algorithm is not support \n',struct_algorithm);
        dag = zeros(r,r);
        return
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% classified test
function [confusion_matrix,correct_rate] = classified_test(dag, classify_node_num,...
    data, classify_algorithm, varargin)

node_flag = 'A';
classify_times = 10;

if ~isempty(varargin{1,1})
args = varargin{1,1}(1,:);
nargs = length(args);
if length(args) > 0
    if ischar(args{1})
        for i=1:2:nargs
            switch args{i}
                case 'node_flag'
                    node_flag = args{i+1};
                case 'classify_times'
                    classify_times = args{i+1};
                case 'params'
                    if isempty(args{i+1}), params = cell(1,n);
                    else params = args{i+1};
                    end
            end
        end
    end
end
end

% set a
ascii_a = 64;

[r,c] = size(data);

confusion_matrix = cell(0,classify_times);
rate = zeros(0,0);
correct_rate =zeros(0,0);

from_classify_node_line = 0;
to_classify_node_line = 0;
for i = 1:r
    if classify_node_num ~= i
        if(dag(classify_node_num,i) == 1)
            from_classify_node_line = from_classify_node_line + 1;
        end
        if(dag(i,classify_node_num) == 1)
            to_classify_node_line = to_classify_node_line + 1;
        end
    end
end

if ~(from_classify_node_line ~= 0 || to_classify_node_line == 0)
    fprintf('dag is not a classifier \n');
    return;
end

% get class_kind
class_kind_map = num2cell(unique(data(:,classify_node_num)));
class_kind = num2cell(unique(data(:,classify_node_num)));
for i = 1:length(class_kind_map)
    class_kind_map{i,2} = num2str(char(ascii_a+i));
    class_kind{i} = num2str(char(ascii_a+i));
end

transposition_data = data';

for i=1:r
    node_data{i} = transposition_data(:,i);
    eval([node_flag,num2str(i), '= transposition_data(:,i);']);
    eval(['stu.',node_flag,num2str(i),'= transposition_data(:,i);']);
    if i == classify_node_num
        classifier_node_data = num2cell(char(transposition_data(:,i) + ascii_a));
        eval([node_flag,num2str(i), '= classifier_node_data;']);
        node_data{i} = classifier_node_data;
        for j = 1:c
            for k = 1:length(class_kind_map)
                if transposition_data(j,i) == cell2num(class_kind_map(k,1))
                    transposition_data(j,i) = cell2num(class_kind_map(k,2));
                end
            end
        end
    end
end

classify_node_num_in_node_data = node_data(1,classify_node_num);
goal = cell2num(classify_node_num_in_node_data);
idx = ismember(cell2num(classify_node_num_in_node_data), class_kind);

stu = structfun(@(x)x(idx,:), stu, 'UniformOutput',false);
numInst = sum(idx);

[goal,~,gnB] = grp2idx(cellstr(cell2num(goal)));

nodeNames = fieldnames(stu);
numNodes = numel(nodeNames);
node = [nodeNames num2cell((1:numNodes)')]';
dNodes = cell2num(node(2,:));
node = struct(node{:});


vals = cell(1,numNodes);
vals(dNodes) = cellfun(@(f) unique(stu.(f)), nodeNames(dNodes), 'Uniform',false);
nodeSize = ones(1,numNodes);
nodeSize(dNodes) = cellfun(@numel, vals(dNodes));

% start classifier
bnet = mk_bnet(dag, nodeSize, 'discrete',dNodes, 'names',nodeNames);
for i=1:numel(dNodes)
    name = nodeNames{dNodes(i)};
    bnet.CPD{dNodes(i)} = tabular_CPD(bnet, node.(name), ...
        'prior_type','dirichlet');
end

%# build samples as cellarray
data = num2cell(cell2mat(struct2cell(stu)')');

% start classifier
for k = 1:classify_times
    
    switch classify_algorithm
        
        case 'HOLD_OUT'
            cv = cvpartition(goal, 'HoldOut',1/3);
        case 'CV-5'
            cv = cvpartition(goal, 'kfold',5);
        case 'CV-10'
            cv = cvpartition(goal, 'kfold',10);
        otherwise
            % error classify algorithm
            fprintf('%s classifier algorithm is not support',classify_algorithm);
            confusion_matrix = [];
            correct_rate = [];
            return;
    end
    
    for i = 1:cv.NumTestSets
        trainData= data(:,cv.training(i));
        testData =data(:,cv.test(i));
        testData(classify_node_num,:) = {[]}; %# remove class
        
        %# training
        bnet = learn_params(bnet, trainData);
        
        %# testing
        prob = zeros(nodeSize(dNodes(1,classify_node_num)), sum(cv.test(i)));
        engine = jtree_inf_engine(bnet);         %# Inference engine
        for j = 1:size(testData,2)
            [engine,loglik] = enter_evidence(engine, testData(:,j));
            marg = marginal_nodes(engine, dNodes(1,classify_node_num));
            prob(:,j) = marg.T;
        end
        
        [~,pred] = max(prob);
        actual= goal(cv.test(i))';
        %# confusion matrix
        predInd = full(sparse(1:numel(pred),pred,1));
        actualInd = full(sparse(1:numel(actual),actual,1));
        [C, RATE] = confmat(predInd, actualInd);
        % confusion_matrix
        confusion_matrix{i,k} = C;
        rate(i) = RATE(:,1);
    end
    
    correct_rate(k) = sum(rate)/cv.NumTestSets;
end

correct_rate = correct_rate';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dag = get_optimal_dag_by_mcmc(data,node_sizes,varargin)

score = 'bic';
nsamples = round(length(data)/5);

args = varargin;
nargs = length(args);
if length(args) > 0
    if ischar(args{1})
        for i=1:2:nargs
            switch args{i}
                case 'score'
                    score = args{i+1};
                case 'nsamples'
                    nsamples = args{i+1};
                case 'params'
                    if isempty(args{i+1})
                        params = cell(1,n);
                    else
                        params = args{i+1};
                    end
            end
        end
    end
end

% mcmc algorithm
dags = learn_struct_mcmc(data,node_sizes,'nsamples',nsamples,...
    'scoring_fn',score);

dag_score = cell(1,length(dags));

for i = 1:length(dags)
    pre_dag = dags{1,i};
    dag_score{i} = score_dags(data, node_sizes, {pre_dag},'scoring_fn',score);
end

% get max score dag
[~,num] = max(cell2mat(dag_score));

dag = dags{1,num};

end