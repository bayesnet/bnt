function [mi] = calculate_mutual_information_array(data)
% FUNCTION [MI_ARRAY] = CALCULATE_MUTUAL_INFORMATION_ARRAY(DATA)
% calculates the mutual information between all pairs of variables
% Data must be discrete, and take values 1,2,...,size
% data(i,m) is the node i in the case m.

[num_nodes num_examples] = size(data);

node_sizes = max(data');
for i = 1:num_nodes
  for ic = 1:node_sizes(i) % I CLASS ic
    px(i,ic) = sum(data(i,:)==ic);
    for j = 1:num_nodes    % J CLASS jc
      for jc = 1:node_sizes(j)
        pxy(i,ic,j,jc) = sum( (data(i,:)==ic) & (data(j,:)==jc) );
      end
      mi(i,j) = 0;
    end
  end
end

for i = 1:num_nodes
  for ic = 1:node_sizes(i)
    for j = 1:num_nodes
      for jc = 1:node_sizes(j)
        if( pxy(i,ic,j,jc)~=0 & px(i,ic)~=0 & px(j,jc)~= 0)
          mi(i,j) = mi(i,j) + pxy(i,ic,j,jc)*log2( num_examples*pxy(i,ic,j,jc)/(px(i,ic)*px(j,jc)) )/num_examples; 
        end
      end
    end
  end
end
