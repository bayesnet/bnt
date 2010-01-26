function slice = get_slice_dbn(bnet, state, i, n, j, m, strides, families, ...
				 CPT)
% slice = get_slice(bnet, state, i, n, j, m, strides, families, cpt)
%
% GET_SLICE get one-dimensional slice of the CPT for node X_i^n
% that corresponds to the different values of X_j^m, where all
% other nodes have values given by state.  
% strides is the result of
% calling compute_strides(bnet)
% families is the result of calling compute_families(bnet)
% cpts is the result of calling get_cpts(bnet)
%
% slice is a 1-d array


if (n == 1)

  k = bnet.eclass1(i);
  c = CPT{k};
  
  % Figure out evidence on family
  fam = families{i, 1};
  ev = state(fam, 1);
  
  % Remove evidence on node j
  pos = find(fam == j);
  ev(pos) = 1;
  dim = size(ev, 1);
  
  % Compute initial index and stride
  start_ind = 1+strides(k, 1:dim)*(ev-1);
  stride = strides(k, pos);

  % Compute the slice
  slice = c(start_ind:stride:start_ind+(bnet.node_sizes(j, 1)-1)*stride);
						  
else
  
  k = bnet.eclass2(i);
  c = CPT{k};
  
  fam = families{i, 2};
  ss = length(bnet.intra);
  
  % Divide the family into nodes in this time step and nodes in the
  % previous time step
  this_time_step = fam(find(fam > ss));
  prev_time_step = fam(find(fam <= ss));

  % Normalize the node numbers
  this_time_step = this_time_step - ss;
  
  % Get the evidence
  this_step_ev = state(this_time_step, n);
  prev_step_ev = state(prev_time_step, n-1);
  
  % Remove the evidence for X_j^m
  if (m == n)
    pos = find(this_time_step == j);
    this_step_ev(pos) = 1;
    pos = pos + size(prev_time_step, 2);
  else
    assert (m == n-1);
    pos = find(prev_time_step == j);
    prev_step_ev(pos) = 1;
  end
  
  % Combine the two time steps
  ev = [prev_step_ev; this_step_ev];
  dim = size(ev, 1);


  % Compute starting index and stride
  start_ind = 1 + strides(k, 1:dim)*(ev-1);
  stride = strides(k, pos);
  
  % Compute slice 
  if (m == 1)
    q = 1;
  else
    q = 2;
  end
  slice = c(start_ind:stride:start_ind+(bnet.node_sizes(j, q)-1)*stride);
end



