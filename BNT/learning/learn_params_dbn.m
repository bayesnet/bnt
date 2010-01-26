function bnet = learn_params_dbn(bnet, data)
% LEARN_PARAM_DBN Estimate params of a DBN for a fully observed model
% bnet = learn_params_dbn(bnet, data)
%
% data(i,t) is the value of node i in slice t (can be a cell array)
% We currently assume there is a single time series
%
% We set bnet.CPD{i} to its ML/MAP estimate.
%
% Currently we assume each node in the first 2 slices has its own CPD (no param tying);
% all nodes in slices >2 share their params with slice 2 as usual.

[ss T] = size(data);

% slice 1
for j=1:ss
  if adjustable_CPD(bnet.CPD{j})
    fam = family(bnet.dag,j);
    bnet.CPD{j} = learn_params(bnet.CPD{j}, data(fam,1));
  end
end


% slices 2:T
% data2(:,t) contains [data(:,t-1); data(:,t)].
% Then we extract out the rows corresponding to the parents in the current and previous slice.
data2 = [data(:,1:T-1);
	 data(:,2:T)];
for j=1:ss
  j2 = j+ss;
  if adjustable_CPD(bnet.CPD{j2})
    fam = family(bnet.dag,j2);
    bnet.CPD{j2} = learn_params(bnet.CPD{j2}, data2(fam,:));
  end
end

