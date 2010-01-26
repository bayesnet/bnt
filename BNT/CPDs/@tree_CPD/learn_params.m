function CPD = learn_params(CPD, fam, data, ns, cnodes, varargin)
% LEARN_PARAMS Construct classification/regression tree given complete data
% CPD = learn_params(CPD, fam, data, ns, cnodes)
%
% fam(i) is the node id of the i-th node in the family of nodes, self node is the last one
% data(i,m) is the value of node i in case m (can be cell array).
% ns(i) is the node size for the i-th node in the whold bnet
% cnodes(i) is the node id for the i-th continuous node in the whole bnet
%  
% The following optional arguments can be specified in the form of name/value pairs:
% stop_cases: for early stop (pruning). A node is not split if it has less than k cases. default is 0.
% min_gain: for early stop (pruning). 
%     For discrete output: A node is not split when the gain of best split is less than min_gain. default is 0.  
%     For continuous (cts) outpt: A node is not split when the gain of best split is less than min_gain*score(root) 
%                                 (we denote it cts_min_gain). default is 0.006
% %%%%%%%%%%%%%%%%%%%Struction definition of dtree_CPD.tree%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tree.num_node               the last position in tree.nodes array for adding new nodes,
%                             it is not always same to number of nodes in a tree, because some position in the 
%                             tree.nodes array can be set to unused (e.g. in tree pruning)  
% tree.nodes is the array of nodes in the tree plus some unused nodes.
% tree.nodes(1) is the root for the tree.
%
% Below is the attributes for each node
% tree.nodes(i).used;     % flag this node is used (0 means node not used, it can be removed from tree to save memory)
% tree.nodes(i).is_leaf;  % if 1 means this node is a leaf, if 0 not a leaf.
% tree.nodes(i).children; % children(i) is the node number in tree.nodes array for the i-th child node
% tree.nodes(i).split_id; % the attribute id used to split this node
% tree.nodes(i).split_threshhold; % the threshhold for continuous attribute to split this node
% %%%%%attributes specially for classification tree (discrete output)
% tree.nodes(i).probs     % probs(i) is the prob for i-th value of class node 
%                         % For three output class, the probs = [0.9 0.1 0.0] means the probability of 
%                         % class 1 is 0.9, for class 2 is 0.1, for class 3 is 0.0.
% %%%%%attributes specially for regression tree (continuous output)                          
% tree.nodes(i).mean      % mean output value for this node
% tree.nodes(i).std       % standard deviation for output values in this node
%
% Author: yimin.zhang@intel.com
% Last updated: Jan. 19, 2002

% Want list:
% (1) more efficient for cts attributes: get the values of cts attributes at first (the begining of build_tree function), then doing bi_search in finding threshhold
% (2) pruning classification tree using Pessimistic Error Pruning
% (3) bi_search for strings (used for transform data to BNT format)

global tree %tree must be global so that it can be accessed in recursive slitting function
global cts_min_gain
tree=[]; % clear the tree
tree.num_node=0;
cts_min_gain=0;

stop_cases=0;
min_gain=0;

args = varargin;
nargs = length(args);
if (nargs>0)
  if isstr(args{1})
    for i=1:2:nargs
      switch args{i},
        case 'stop_cases', stop_cases = args{i+1};   
        case 'min_gain', min_gain = args{i+1};
      end
    end
  else
    error(['error in input parameters']);
  end
end

if iscell(data)
  local_data = cell2num(data(fam,:));
else
  local_data = data(fam, :);
end
%counts = compute_counts(local_data, CPD.sizes);
%CPD.CPT = mk_stochastic(counts + CPD.prior); % bug fix 11/5/01
node_types = zeros(1,size(ns,2)); %all nodes are disrete
node_types(cnodes)=1;
%make the data be BNT compliant (values for discrete nodes are from 1-n, here n is the node size)
%trans_data=transform_data(local_data,'tmp.dat',[]); %here no cts nodes

build_dtree (CPD, local_data, ns(fam), node_types(fam),stop_cases,min_gain);
%CPD.tree=copy_tree(tree);
CPD.tree=tree; %copy the tree constructed to CPD


function new_tree = copy_tree(tree)
% copy the tree to new_tree
new_tree.num_node=tree.num_node;
new_tree.root = tree.root;
for i=1:tree.num_node
  new_tree.nodes(i)=tree.nodes(i);
end


function build_dtree (CPD, fam_ev, node_sizes, node_types,stop_cases,min_gain)
global tree
global cts_min_gain

tree.num_node=0; %the current number of nodes in the tree
tree.root=1;

T = 1:size(fam_ev,2) ; %all cases
candidate_attrs = 1:(size(node_sizes,2)-1); %all attributes
node_id=1;  %the root node
lastnode=size(node_sizes,2); %the last element in all nodes is the dependent variable (category node)
num_cat=node_sizes(lastnode);

% get minimum gain for cts output (used in stop splitting)
if (node_types(size(fam_ev,1))==1) %cts output
  N = size(fam_ev,2);
  output_id = size(fam_ev,1);
  cases_T = fam_ev(output_id,:); %get all the output value for cases T
  std_T = std(cases_T);
  avg_y_T = mean(cases_T);
  sqr_T = cases_T - avg_y_T;
  cts_min_gain = min_gain*(sum(sqr_T.*sqr_T)/N);  % min_gain * (R(root) = 1/N * SUM(y-avg_y)^2)
end  

split_dtree (CPD, fam_ev, node_sizes, node_types, stop_cases,min_gain, T, candidate_attrs, num_cat);
  


% pruning method
% (1) Restrictions on minimum node size: A node is not split if it has smaller than k cases.
% (2) Threshholds on impurity: a threshhold is imposed on the splitting test score. Threshhold can be 
% imposed on local goodness measure (the gain_ratio of a node) or global goodness.
% (3) Mininum Error Pruning (MEP), (no need pruning set)
%     Prune if static error<=backed-up error
%      Static error at node v: e(v) = (Nc + 1)/(N+k) (laplace estimate, prior for each class equal) 
%        here N is # of all examples, Nc is # of majority class examples, k is number of classes 
%      Backed-up error at node v: (Ti is the i-th subtree root)
%         E(T) = Sum_1_to_n(pi*e(Ti))
% (4) Pessimistic Error Pruning (PEP), used in Quilan C4.5 (no need pruning set, efficient because of pruning top-down)
%       Probability of error (apparent error rate)
%           q = (N-Nc+0.5)/N
%         where N=#examples, Nc=#examples in majority class
%     Error of a node v (if pruned)  q(v)= (Nv- Nc,v + 0.5)/Nv
%     Error of a subtree   q(T)= Sum_of_l_leaves(Nl - Nc,l + 0.5)/Sum_of_l_leaves(Nl)
%     Prune if q(v)<=q(T)
% 
% Implementation statuts:
% (1)(2) has been implemented as the input parameters of learn_params.
% (4) is implemented in this function
function pruning(fam_ev,node_sizes,node_types)
% PRUNING prune the constructed tree using PEP
% pruning(fam_ev,node_sizes,node_types)
%
% fam_ev(i,j)  is the value of attribute i in j-th training cases (for whole tree), the last row is for the class label (self_ev)
% node_sizes(i) is the node size for the i-th node in the family
% node_types(i) is the node type for the i-th node in the family, 0 for disrete node, 1 for continous node
% the global parameter 'tree' is for storing the input tree and the pruned tree


function split_T = split_cases(fam_ev,node_sizes,node_types,T,node_i, threshhold)
% SPLIT_CASES split the cases T according to values of node_i in the family
% split_T = split_cases(fam_ev,node_sizes,node_types,T,node_i)
%
% fam_ev(i,j)  is the value of attribute i in j-th training cases (for whole tree), the last row is for the class label (self_ev)
% node_sizes(i) is the node size for the i-th node in the family
% node_types(i) is the node type for the i-th node in the family, 0 for disrete node, 1 for continous node
% node_i is the attribute we need to split

if (node_types(node_i)==0) %discrete attribute
  %init the subsets of T
  split_T = cell(1,node_sizes(node_i)); %T will be separated into |node_size of i| subsets according to different values of node i
  for i=1:node_sizes(node_i)   % here we assume that the value of an attribute is 1:node_size
    split_T{i}=zeros(1,0);
  end

  size_t = size(T,2);
  for i=1:size_t
    case_id = T(i);
    %put this case into one subset of split_T according to its value for node_i
    value = fam_ev(node_i,case_id); 
    pos = size(split_T{value},2)+1;
    split_T{value}(pos)=case_id;  % here assumes the value of an attribute is 1:node_size 
  end
else %continuous attribute
  %init the subsets of T
  split_T = cell(1,2); %T will be separated into 2 subsets (<=threshhold) (>threshhold)
  for i=1:2   
    split_T{i}=zeros(1,0);
  end

  size_t = size(T,2);
  for i=1:size_t
    case_id = T(i);
    %put this case into one subset of split_T according to its value for node_i
    value = fam_ev(node_i,case_id); 
    subset_num=1;
    if (value>threshhold)
      subset_num=2;
    end  
    pos = size(split_T{subset_num},2)+1;
    split_T{subset_num}(pos)=case_id;  
  end
end


  
function new_node = split_dtree (CPD, fam_ev, node_sizes, node_types, stop_cases, min_gain, T, candidate_attrs, num_cat)
% SPLIT_TREE Split the tree at node node_id with cases T (actually it is just indexes to family evidences).
% new_node = split_dtree (fam_ev, node_sizes, node_types, T, node_id, num_cat, method)
%
% fam_ev(i,j)  is the value of attribute i in j-th training cases (for whole tree), the last row is for the class label (self_ev)
% node_sizes{i} is the node size for the i-th node in the family
% node_types{i} is the node type for the i-th node in the family, 0 for disrete node, 1 for continous node
% stop_cases is the threshold of number of cases to stop slitting
% min_gain is the minimum gain need to split a node
% T(i) is the index of i-th cases in current decision tree node, we need split it further
% candidate_attrs(i) the node id for the i-th attribute that still need to be considered as split attribute 
%%%%% node_id is the index of current node considered for a split
% num_cat is the number of output categories for the decision tree
% output:
% new_node is the new node created
global tree
global cts_min_gain

size_fam = size(fam_ev,1);            %number of family size
output_type = node_types(size_fam);   %the type of output for the tree (0 is discrete, 1 is continuous)
size_attrs = size(candidate_attrs,2); %number of candidate attributes
size_t = size(T,2);                   %number of training cases in this tree node

%(1)computeFrequenceyForEachClass(T)
if (output_type==0) %discrete output
  class_freqs = zeros(1,num_cat);
  for i=1:size_t
    case_id = T(i);
    case_class = fam_ev(size_fam,case_id); %get the class label for this case
    class_freqs(case_class)=class_freqs(case_class)+1;
  end
else  %cts output
  N = size(fam_ev,2);
  cases_T = fam_ev(size(fam_ev,1),T); %get the output value for cases T
  std_T = std(cases_T);
end

%(2) if OneClass (for discrete output) or same output value (for cts output) or Class With #examples < stop_cases
%         return a leaf;
%    create a decision node N;

% get majority class in this node
if (output_type == 0)
  top1_class = 0;       %the class with the largest number of cases
  top1_class_cases = 0; %the number of cases in top1_class
  [top1_class_cases,top1_class]=max(class_freqs);
end
  
if (size_t==0)     %impossble
  new_node=-1;
  fprintf('Fatal error: please contact the author. \n');
  return;
end

% stop splitting if needed
  %for discrete output: one class 
  %for cts output, all output value in cases are same
  %cases too little
if ( (output_type==0 & top1_class_cases == size_t) | (output_type==1 & std_T == 0) | (size_t < stop_cases))             
  %create one new leaf node
  tree.num_node=tree.num_node+1;
  tree.nodes(tree.num_node).used=1; %flag this node is used (0 means node not used, it will be removed from tree at last to save memory)
  tree.nodes(tree.num_node).is_leaf=1;
  tree.nodes(tree.num_node).children=[];
  tree.nodes(tree.num_node).split_id=0;  %the attribute(parent) id to split this tree node
  tree.nodes(tree.num_node).split_threshhold=0;  
  if (output_type==0)
    tree.nodes(tree.num_node).probs=class_freqs/size_t; %the prob for each value of class node 

    %  tree.nodes(tree.num_node).probs=zeros(1,num_cat); %the prob for each value of class node 
    %  tree.nodes(tree.num_node).probs(top1_class)=1; %use the majority class of parent node, like for binary class, 
                                                   %and majority is class 2, then the CPT is [0 1]
                                                   %we may need to use prior to do smoothing, to get [0.001 0.999]
    tree.nodes(tree.num_node).error.self_error=1-top1_class_cases/size_t; %the classfication error in this tree node when use default class
    tree.nodes(tree.num_node).error.all_error=1-top1_class_cases/size_t;  %no total classfication error in this tree node and its subtree
    tree.nodes(tree.num_node).error.all_error_num=size_t - top1_class_cases;
    fprintf('Create leaf node(onecla) %d. Class %d Cases %d Error %d \n',tree.num_node, top1_class, size_t, size_t - top1_class_cases );
  else
    avg_y_T = mean(cases_T);
    tree.nodes(tree.num_node).mean = avg_y_T; 
    tree.nodes(tree.num_node).std = std_T;
    fprintf('Create leaf node(samevalue) %d. Mean %8.4f Std %8.4f Cases %d \n',tree.num_node, avg_y_T, std_T, size_t);
  end  
  new_node = tree.num_node;
  return;
end
    
%create one new node
tree.num_node=tree.num_node+1;
tree.nodes(tree.num_node).used=1; %flag this node is used (0 means node not used, it will be removed from tree at last to save memory)
tree.nodes(tree.num_node).is_leaf=1;
tree.nodes(tree.num_node).children=[];
tree.nodes(tree.num_node).split_id=0;
tree.nodes(tree.num_node).split_threshhold=0;  
if (output_type==0)
  tree.nodes(tree.num_node).error.self_error=1-top1_class_cases/size_t; 
  tree.nodes(tree.num_node).error.all_error=0;
  tree.nodes(tree.num_node).error.all_error_num=0;
else
  avg_y_T = mean(cases_T);
  tree.nodes(tree.num_node).mean = avg_y_T; 
  tree.nodes(tree.num_node).std = std_T;
end
new_node = tree.num_node;

%Stop splitting if no attributes left in this node
if (size_attrs==0) 
  if (output_type==0)
    tree.nodes(tree.num_node).probs=class_freqs/size_t; %the prob for each value of class node 
    tree.nodes(tree.num_node).error.all_error=1-top1_class_cases/size_t;  
    tree.nodes(tree.num_node).error.all_error_num=size_t - top1_class_cases;
    fprintf('Create leaf node(noattr) %d. Class %d Cases %d Error %d \n',tree.num_node, top1_class, size_t, size_t - top1_class_cases );
  else
    fprintf('Create leaf node(noattr) %d. Mean %8.4f Std %8.4f Cases %d \n',tree.num_node, avg_y_T, std_T, size_t);
  end
  return;
end
      
  
%(3) for each attribute A
%        ComputeGain(A);
max_gain=0;  %the max gain score (for discrete information gain or gain ration, for cts node the R(T))
best_attr=0;  %the attribute with the max_gain
best_split = []; %the split of T according to the value of best_attr
cur_best_threshhold = 0; %the threshhold for split continuous attribute
best_threshhold=0;

% compute Info(T) (for discrete output)
if (output_type == 0)
  class_split_T = split_cases(fam_ev,node_sizes,node_types,T,size(fam_ev,1),0); %split cases according to class
  info_T = compute_info (fam_ev, T, class_split_T);
else % compute R(T) (for cts output)
%  N = size(fam_ev,2);
%  cases_T = fam_ev(size(fam_ev,1),T); %get the output value for cases T
%  std_T = std(cases_T);
%  avg_y_T = mean(cases_T);
  sqr_T = cases_T - avg_y_T;
  R_T = sum(sqr_T.*sqr_T)/N;  % get R(T) = 1/N * SUM(y-avg_y)^2
  info_T = R_T;
end

for i=1:(size_fam-1)
  if (myismember(i,candidate_attrs))  %if this attribute still in the candidate attribute set
    if (node_types(i)==0) %discrete attibute
      split_T = split_cases(fam_ev,node_sizes,node_types,T,i,0); %split cases according to value of attribute i
      % For cts output, we compute the least square gain.
      % For discrete output, we compute gain ratio
      cur_gain = compute_gain(fam_ev,node_sizes,node_types,T,info_T,i,split_T,0,output_type); %gain ratio
    else %cts attribute
      %get the values of this attribute
      ev = fam_ev(:,T);
      values = ev(i,:);
      sort_v = sort(values); 
        %remove the duplicate values in sort_v
      v_set = unique(sort_v);  
      best_gain = 0;
      best_threshhold = 0;
      best_split1 = [];
      
      %find the best split for this cts attribute
      % see "Quilan 96: Improved Use of Continuous Attributes in C4.5"
      for j=1:(size(v_set,2)-1)
        mid_v = (v_set(j)+v_set(j+1))/2; 
        split_T = split_cases(fam_ev,node_sizes,node_types,T,i,mid_v); %split cases according to value of attribute i (<=mid_v)
        % For cts output, we compute the least square gain.
        % For discrete output, we use Quilan 96: use information gain instead of gain ratio to select threshhold
        cur_gain = compute_gain(fam_ev,node_sizes,node_types,T,info_T,i,split_T,1,output_type); 
        %if (i==6)
        %  fprintf('gain %8.5f threshhold %6.3f spliting %d\n', cur_gain, mid_v, size(split_T{1},2));
        %end

        if (best_gain < cur_gain)
          best_gain = cur_gain;
          best_threshhold = mid_v;
          %best_split1 = split_T;     %here we need to copy array, not good!!! (maybe we can compute after we get best_attr
        end
      end
      %recalculate the gain_ratio of the best_threshhold
      split_T = split_cases(fam_ev,node_sizes,node_types,T,i,best_threshhold);
      best_gain = compute_gain(fam_ev,node_sizes,node_types,T,info_T,i,split_T,0,output_type); %gain_ratio
      if (output_type==0) %for discrete output
        cur_gain = best_gain-log2(size(v_set,2)-1)/size_t; % Quilan 96: use the gain_ratio-log2(N-1)/|D| as the gain of this attr
      else                %for cts output
        cur_gain = best_gain;
      end
    end
    
    if (max_gain < cur_gain)
      max_gain = cur_gain;
      best_attr = i;
      cur_best_threshhold=best_threshhold;  %save the threshhold
      %best_split = split_T;        %here we need to copy array, not good!!! So we will recalculate in below line 313
    end
  end
end

% stop splitting if gain is too small
if (max_gain==0 | (output_type==0 & max_gain < min_gain) | (output_type==1 & max_gain < cts_min_gain)) 
  if (output_type==0)
    tree.nodes(tree.num_node).probs=class_freqs/size_t; %the prob for each value of class node 
    tree.nodes(tree.num_node).error.all_error=1-top1_class_cases/size_t;  
    tree.nodes(tree.num_node).error.all_error_num=size_t - top1_class_cases;
    fprintf('Create leaf node(nogain) %d. Class %d Cases %d Error %d \n',tree.num_node, top1_class, size_t, size_t - top1_class_cases );
  else
    fprintf('Create leaf node(nogain) %d. Mean %8.4f Std %8.4f Cases %d \n',tree.num_node, avg_y_T, std_T, size_t);
  end
  return;
end

%get the split of cases according to the best split attribute
if (node_types(best_attr)==0) %discrete attibute
  best_split = split_cases(fam_ev,node_sizes,node_types,T,best_attr,0);  
else  
  best_split = split_cases(fam_ev,node_sizes,node_types,T,best_attr,cur_best_threshhold);
end
  
%(4) best_attr = AttributeWithBestGain;
%(5) if best_attr is continuous             ???? why need this? maybe the value in the decision tree must appeared in data
%       find threshhold in all cases that <= max_V
%    change the split of T
tree.nodes(tree.num_node).split_id=best_attr;
tree.nodes(tree.num_node).split_threshhold=cur_best_threshhold; %for cts attribute only

%note: below threshhold rejust is linera search, so it is slow. A better method is described in paper "Efficient C4.5"
%if (output_type==0)
if (node_types(best_attr)==1)  %is a continuous attribute
  %find the value that approximate best_threshhold from below (the largest that <= best_threshhold)
  best_value=0;
  for i=1:size(fam_ev,2)  %note: need to search in all cases for all tree, not just in cases for this node
    val = fam_ev(best_attr,i);
    if (val <= cur_best_threshhold & val > best_value) %val is more clear to best_threshhold
      best_value=val;
    end
  end
  tree.nodes(tree.num_node).split_threshhold=best_value; %for cts attribute only
end
%end
  
if (output_type == 0)
  fprintf('Create node %d split at %d gain %8.4f Th %d. Class %d Cases %d Error %d \n',tree.num_node, best_attr, max_gain, tree.nodes(tree.num_node).split_threshhold, top1_class, size_t, size_t - top1_class_cases );
else
  fprintf('Create node %d split at %d gain %8.4f Th %d. Mean %8.4f Cases %d\n',tree.num_node, best_attr, max_gain, tree.nodes(tree.num_node).split_threshhold, avg_y_T, size_t );
end
  
%(6) Foreach T' in the split_T
%        if T' is Empty
%            Child of node_id is a leaf
%        else
%            Child of node_id = split_tree (T')
tree.nodes(new_node).is_leaf=0; %because this node will be split, it is not leaf now
for i=1:size(best_split,2)
  if (size(best_split{i},2)==0) %T(i) is empty
    %create one new leaf node
    tree.num_node=tree.num_node+1;
    tree.nodes(tree.num_node).used=1; %flag this node is used (0 means node not used, it will be removed from tree at last to save memory)
    tree.nodes(tree.num_node).is_leaf=1;
    tree.nodes(tree.num_node).children=[];
    tree.nodes(tree.num_node).split_id=0;
    tree.nodes(tree.num_node).split_threshhold=0;  
    if (output_type == 0)
      tree.nodes(tree.num_node).probs=zeros(1,num_cat); %the prob for each value of class node 
      tree.nodes(tree.num_node).probs(top1_class)=1; %use the majority class of parent node, like for binary class, 
                                                   %and majority is class 2, then the CPT is [0 1]
                                                   %we may need to use prior to do smoothing, to get [0.001 0.999]
      tree.nodes(tree.num_node).error.self_error=0; 
      tree.nodes(tree.num_node).error.all_error=0;  
      tree.nodes(tree.num_node).error.all_error_num=0;
    else
      tree.nodes(tree.num_node).mean = avg_y_T; %just use parent node's mean value
      tree.nodes(tree.num_node).std = std_T;
    end
    %add the new leaf node to parents
    num_children=size(tree.nodes(new_node).children,2);
    tree.nodes(new_node).children(num_children+1)=tree.num_node;
    if (output_type==0)
      fprintf('Create leaf node(nullset) %d. %d-th child of Father %d Class %d\n',tree.num_node, i, new_node, top1_class );
    else
      fprintf('Create leaf node(nullset) %d. %d-th child of Father %d \n',tree.num_node, i, new_node );
    end

  else
    if (node_types(best_attr)==0)  % if attr is discrete, it should be removed from the candidate set  
      new_candidate_attrs = mysetdiff(candidate_attrs,[best_attr]);
    else
      new_candidate_attrs = candidate_attrs;
    end
    new_sub_node = split_dtree (CPD, fam_ev, node_sizes, node_types, stop_cases, min_gain, best_split{i}, new_candidate_attrs, num_cat);  
    %tree.nodes(parent_id).error.all_error += tree.nodes(new_sub_node).error.all_error;
    fprintf('Add subtree node %d to %d. #nodes %d\n',new_sub_node,new_node, tree.num_node );

%   tree.nodes(new_node).error.all_error_num = tree.nodes(new_node).error.all_error_num + tree.nodes(new_sub_node).error.all_error_num;
    %add the new leaf node to parents
    num_children=size(tree.nodes(new_node).children,2);
    tree.nodes(new_node).children(num_children+1)=new_sub_node;
  end
end   
  
%(7) Compute errors of N; for doing pruning
%    get the total error for the subtree
if (output_type==0)
  tree.nodes(new_node).error.all_error=tree.nodes(new_node).error.all_error_num/size_t;
end
%doing pruning, but doing here is not so efficient, because it is bottom up.
%if tree.nodes()
%after doing pruning, need to update the all_error to self_error

%(8) Return N
  



%(1) For discrete output, we use GainRatio defined as below
%  			         Gain(X,T)
% 	GainRatio(X,T) = ----------
% 			         SplitInfo(X,T)
%   where
%   Gain(X,T) = Info(T) - Info(X,T)
%    				                       |Ti|
% 	Info(X,T) = Sum for i from 1 to n of ( ---- * Info(Ti))
%                                          |T|
 			 
%   SplitInfo(D,T) is the information due to the split of T on the basis
%    of the value of the categorical attribute D. Thus SplitInfo(D,T) is
%  		 I(|T1|/|T|, |T2|/|T|, .., |Tm|/|T|)
%    where {T1, T2, .. Tm} is the partition of T induced by the value of D.

%   Definition of Info(Ti)
%     If a set T of records is partitioned into disjoint exhaustive classes C1, C2, .., Ck on the basis of the 
%     value of the categorical attribute, then the information needed to identify the class of an element of T 
%     is Info(T) = I(P), where P is the probability distribution of the partition (C1, C2, .., Ck): 
%     	P = (|C1|/|T|, |C2|/|T|, ..., |Ck|/|T|)
%     Here I(P) is defined as
%       I(P) = -(p1*log(p1) + p2*log(p2) + .. + pn*log(pn))
% 
%(2) For continuous output (regression tree), we use least squares score (adapted from Leo Breiman's book "Classification and regression trees", page 231
%    The original support only binary split, we further extend it to permit multiple-child split
%                                        
%     Delta_R = R(T) - Sum for all childe nodes Ti (R(Ti))
%     Where R(Ti)= 1/N * Sum for all cases i in node Ti ((yi - avg_y(Ti))^2)
%     here N is the number of all training cases for construct the regression tree
%          avg_y(Ti) is the average value for output variable for the cases in node Ti

function gain_score = compute_gain (fam_ev, node_sizes, node_types, T, info_T, attr_id, split_T, score_type, output_type)
% COMPUTE_GAIN Compute the score for the split of cases T using attribute attr_id
% gain_score = compute_gain (fam_ev, T, attr_id, node_size, method)
%
% fam_ev(i,j)  is the value of attribute i in j-th training cases, the last row is for the class label (self_ev)
% T(i) is the index of i-th cases in current decision tree node, we need split it further
% attr_id is the index of current node considered for a split
% split_T{i} is the i_th subset in partition of cases T according to the value of attribute attr_id
% score_type if 0, is gain ratio, 1 is information gain (only apply to discrete output)
% node_size(i) the node size of i-th node in the family
% output_type: 0 means discrete output, 1 means continuous output.
gain_score=0;
% ***********for DISCRETE output*******************************************************
if (output_type == 0)
  % compute Info(T)
  total_cnt = size(T,2);
  if (total_cnt==0)
    return;
  end;
  %class_split_T = split_cases(fam_ev,node_sizes,node_types,T,size(fam_ev,1),0); %split cases according to class
  %info_T = compute_info (fam_ev, T, class_split_T);

  % compute Info(X,T)
  num_class = size(split_T,2); 
  subset_sizes = zeros(1,num_class);
  info_ti = zeros(1,num_class);
  for i=1:num_class
    subset_sizes(i)=size(split_T{i},2);
    if (subset_sizes(i)~=0)
      class_split_Ti = split_cases(fam_ev,node_sizes,node_types,split_T{i},size(fam_ev,1),0); %split cases according to class
      info_ti(i) = compute_info(fam_ev, split_T{i}, class_split_Ti);
    end
  end    
  ti_ratios = subset_sizes/total_cnt;  %get the |Ti|/|T|
  info_X_T = sum(ti_ratios.*info_ti);

  %get Gain(X,T)
  gain_X_T = info_T - info_X_T;

  if (score_type == 1) %information gain
    gain_score=gain_X_T;
    return;
  end
  %compute the SplitInfo(X,T)   //is this also for cts attr, only split into two subsets
  splitinfo_T = compute_info (fam_ev, T, split_T);
  if (splitinfo_T~=0)
    gain_score = gain_X_T/splitinfo_T;
  end

% ************for continuous output**************************************************
else 
  N = size(fam_ev,2);

  % compute R(Ti)
  num_class = size(split_T,2); 
  R_Ti = zeros(1,num_class);
  for i=1:num_class
    if (size(split_T{i},2)~=0)
      cases_T = fam_ev(size(fam_ev,1),split_T{i});
      avg_y_T = mean(cases_T);
      sqr_T = cases_T - avg_y_T;
      R_Ti(i) = sum(sqr_T.*sqr_T)/N;  % get R(Ti) = 1/N * SUM(y-avg_y)^2
    end
  end
  %delta_R = R(T) - SUM(R(Ti))
  gain_score = info_T - sum(R_Ti);

end


%   Definition of Info(Ti)
%     If a set T of records is partitioned into disjoint exhaustive classes C1, C2, .., Ck on the basis of the 
%     value of the categorical attribute, then the information needed to identify the class of an element of T 
%     is Info(T) = I(P), where P is the probability distribution of the partition (C1, C2, .., Ck): 
%     	P = (|C1|/|T|, |C2|/|T|, ..., |Ck|/|T|)
%     Here I(P) is defined as
%       I(P) = -(p1*log(p1) + p2*log(p2) + .. + pn*log(pn))
function info = compute_info (fam_ev, T, split_T)
% COMPUTE_INFO compute the information for the split of T into split_T
% info = compute_info (fam_ev, T, split_T)

total_cnt = size(T,2);
num_class = size(split_T,2);
subset_sizes = zeros(1,num_class);
probs = zeros(1,num_class);
log_probs = zeros(1,num_class);
for i=1:num_class
  subset_sizes(i)=size(split_T{i},2);
end    

probs = subset_sizes/total_cnt;
%log_probs = log2(probs);  % if probs(i)=0, the log2(probs(i)) will be Inf
for i=1:size(probs,2)
  if (probs(i)~=0)
    log_probs(i)=log2(probs(i));
  end
end

info = sum(-(probs.*log_probs));

