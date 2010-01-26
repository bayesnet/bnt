% Here the training data is adapted from UCI ML repository, 'zoo' data

dtreeCPD=tree_CPD;

% load data
fname = fullfile(BNT_HOME, 'examples', 'static', 'uci_data', 'zoo', 'zoo1.data')
data=load(fname);
data=data';

data=transform_data_into_bnt_format(data, []);

% learn decision tree from data 
ns=2*ones(1,17);
ns(13)=6;
ns(17)=7;
dtreeCPD1=learn_params(dtreeCPD,1:17,data,ns,[],'stop_cases',5); % a node with less than 5 cases will not be splitted

% evaluate on data
[score,outputs]=evaluate_tree_performance(dtreeCPD1,1:17,data,ns,[]);
fprintf('Accuracy in old training data %6.3f\n',score);

