function [bnt_data, old_values] = transform_data_into_bnt_format(data,cnodes)
% TRANSFORM_DATA_TO_BNT_FORMAT Ensures discrete variables have values 1,2,..,k
% e.g., if the values of a discrete are [0 1 6], they must be mapped to [1 2 3]
%
% data(i,j) is the value for i-th node in j-th case.
% bnt_data(i,j) is the new value.
% old_values{i} are the original values for node i.
% cnodes is the list of all continous nodes, e.g. [3 5] means the 3rd and 5th node is continuous
%
% Author: yimin.zhang@intel.com
% Last updated: Jan. 22, 2002 by Kevin Murphy.

num_nodes=size(data,1);
num_cases=size(data,2);
old_values=cell(1,num_nodes);

for i=1:num_nodes
  if (myismember(i,cnodes)==1)  %cts nodes no need to be transformed 
    %just copy the data
    bnt_data(i,:)=data(i,:);
    continue;
  end
  values = data(i,:);
  sort_v = sort(values); 
  %remove the duplicate values in sort_v
  v_set = unique(sort_v);  
  
  %transform the values
  for j=1:size(values,2)
    index = binary_search(v_set,values(j));
    if (index==-1)
      fprintf('value not found in tranforming data to bnt format.\n');   
      return;
    end
    bnt_data(i,j)=index;
  end
  old_values{i}=v_set;
end


%%%%%%%%%%%%

function index=binary_search(vector, value)
% BI_SEARCH do binary search for value in the vector
% Author: yimin.zhang@intel.com
% Last updated: Jan. 19, 2002

begin_index=1;
end_index=size(vector,2); 
index=-1;
while (begin_index<=end_index)
  mid=floor((begin_index+end_index)/2);
  if (isstr(vector(mid)))
    % need to write a strcmp to return three result (< = >)
  else
    if (value==vector(mid))
      index=mid;
      return;
    elseif (value>vector(mid))
      begin_index=mid+1;    
    else
      end_index=mid-1;
    end
  end
end
return;
