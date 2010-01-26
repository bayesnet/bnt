function [score,outputs] = evaluate(CPD, fam, data, ns, cnodes)
% Evaluate evaluate the performance of the classification/regression tree on given complete data
% score = evaluate(CPD, fam, data, ns, cnodes)
%
% fam(i) is the node id of the i-th node in the family of nodes, self node is the last one
% data(i,m) is the value of node i in case m (can be cell array).
% ns(i) is the node size for the i-th node in the whold bnet
% cnodes(i) is the node id for the i-th continuous node in the whole bnet
%  
% Output
% score is the classification accuracy (for classification) 
%          or mean square deviation (for regression)
%            here for every case we use the mean value at the tree leaf node as its predicted value
% outputs(i) is the predicted output value for case i
%
% Author: yimin.zhang@intel.com
% Last updated: Jan. 19, 2002


if iscell(data)
  local_data = cell2num(data(fam,:));
else
  local_data = data(fam, :);
end

%get local node sizes and node types
node_sizes = ns(fam);
node_types = zeros(1,size(ns,2)); %all nodes are disrete
node_types(cnodes)=1;
node_types=node_types(fam);

fam_size=size(fam,2);
output_type = node_types(fam_size);

num_cases=size(local_data,2);
total_error=0;

outputs=zeros(1,num_cases);
for i=1:num_cases
  %class one case using the tree
  cur_node=CPD.tree.root;  % at the root node of the tree
  while (1)
    if (CPD.tree.nodes(cur_node).is_leaf==1)
      if (output_type==0) %output is discrete
        %use the class with max probability as the output  
        [maxvalue,class_id]=max(CPD.tree.nodes(cur_node).probs);
        outputs(i)=class_id;
        if (class_id~=local_data(fam_size,i))
          total_error=total_error+1;
        end
      else   %output is continuous
        %use the mean as the value
        outputs(i)=CPD.tree.nodes(cur_node).mean;
        cur_deviation = CPD.tree.nodes(cur_node).mean-local_data(fam_size,i);
        total_error=total_error+cur_deviation*cur_deviation;
      end
      break;
    end
    cur_attr = CPD.tree.nodes(cur_node).split_id; 
    attr_val = local_data(cur_attr,i);
    if (node_types(cur_attr)==0)  %discrete attribute
        % goto the attr_val -th child
        cur_node = CPD.tree.nodes(cur_node).children(attr_val);
    else
        if (attr_val <= CPD.tree.nodes(cur_node).split_threshhold)
          cur_node = CPD.tree.nodes(cur_node).children(1);
        else
          cur_node = CPD.tree.nodes(cur_node).children(2);  
        end
    end
    if (cur_node > CPD.tree.num_node)
      fprintf('Fatal error: Tree structure corrupted.\n');
      return;
    end
  end
  %update the classification error number
end
if (output_type==0)
  score=1-total_error/num_cases;
else
  score=total_error/num_cases;
end
