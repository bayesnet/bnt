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
  e = bnet.equiv_class(j);
  if adjustable_CPD(bnet.CPD{e})
    fam = family(bnet.dag,e);
    bnet.CPD{j} = learn_params(bnet.CPD{e}, fam, data, bnet.node_sizes, bnet.cnodes); % F.Denk 01.07.2017 - learn_params input parameters in previous script did not match the required ones
  end
end


% slices 2:T
% data2(:,t) contains [data(:,t-1); data(:,t)].
% Then we extract out the rows corresponding to the parents in the current and previous slice.
data2 = [data(:,1:T-1);
	 data(:,2:T)];
for j=1:ss
  j2 = j+ss;
  e = bnet.equiv_class(j2);
  if adjustable_CPD(bnet.CPD{e})
    fam = family(bnet.dag,j2); %% F.Denk 01.07.2017 family of nodes in slice 2 are same as for nodes in slice 1 except j2 = j + ss -- equivalence class ensures correct use of parameters
    bnet.CPD{e} = learn_params(bnet.CPD{e}, fam, data2, bnet.node_sizes, bnet.cnodes); % F.Denk 01.07.2017 - learn_params input parameters in previous script did not match the required ones
  end
end

