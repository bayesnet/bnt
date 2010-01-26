function [mpe, ll] = calc_mpe_bucket(bnet, new_evidence, max_over)
%
% PURPOSE:
%       CALC_MPE Computes the most probable explanation to the network nodes
%       given the evidence.
%       
%       [mpe, ll] = calc_mpe(engine, new_evidence, max_over)
%
% INPUT:
%       bnet  - the bayesian network
%       new_evidence - optional, if specified - evidence to be incorporated [cell(1,n)]
%       max_over - optional, if specified determines the variable elimination order [1:n]
%
% OUTPUT:
%       mpe - the MPE assignmet for the net variables (or [] if no satisfying assignment)
%       ll - log assignment probability.
%
% Notes:
% 1. Adapted from '@var_elim_inf_engine\marginal_nodes' for MPE by Ron Zohar, 8/7/01
% 2. Only discrete potentials are supported at this time.
% 3. Complexity: O(nw*) where n is the number of nodes and w* is the induced tree width.
% 4. Implementation based on:
%  - R. Dechter, "Bucket Elimination: A Unifying Framework for Probabilistic Inference", 
%                 UA1 96, pp. 211-219.


ns = bnet.node_sizes;
n = length(bnet.dag);
evidence = cell(1,n);
if (nargin<2)
    new_evidence = evidence;
end

onodes = find(~isemptycell(new_evidence));  % observed nodes
hnodes = find(isemptycell(new_evidence));  % hidden nodes
pot_type = determine_pot_type(bnet, onodes);

if pot_type ~= 'd'
  error('only disrete potentials supported at this time')    
end

for i=1:n
  fam = family(bnet.dag, i);
  CPT{i} = convert_to_pot(bnet.CPD{bnet.equiv_class(i)}, pot_type, fam(:), evidence);        
end 

% handle observed nodes: set impossible cases' probability to zero
% rather than prun matrix (this makes backtracking easier)

for ii=onodes
  lIdx = 1:ns(ii);
  lIdx = setdiff(lIdx, new_evidence{ii});
  
  sCPT=struct(CPT{ii});  % violate object privacy
  
  sargs = '';
  for jj=1:(length(sCPT.domain)-1)
    sargs = [sargs, ':,']; 
  end        
  for jj=lIdx
    eval(['sCPT.T(', sargs, num2str(jj), ')=0;']);
  end
  CPT{ii}=dpot(sCPT.domain, sCPT.sizes, sCPT.T);        
end

B = cell(1,n); 
for b=1:n
  B{b} = mk_initial_pot(pot_type, [], [], [], []);
end

if (nargin<3)
  max_over = (1:n);
end   
order = max_over; % no attempt to optimize this


% Initialize the buckets with the CPDs assigned to them
for i=1:n
  b = bucket_num(domain_pot(CPT{i}), order);
  B{b} = multiply_pots(B{b}, CPT{i});
end

% Do backward phase
max_over = max_over(length(max_over):-1:1); % reverse
for i=max_over(1:end-1)        
  % max-ing over variable i which occurs in bucket j
  j = bucket_num(i, order);
  rest = mysetdiff(domain_pot(B{j}), i);
  %temp = marginalize_pot_max(B{j}, rest);
  temp = marginalize_pot(B{j}, rest, 1);
  b = bucket_num(domain_pot(temp), order);
  %        fprintf('maxing over bucket %d (var %d), putting result into bucket %d\n', j, i, b);
  sB=struct(B{b});  % violate object privacy
  if ~isempty(sB.domain)
    B{b} = multiply_pots(B{b}, temp);
  else
    B{b} = temp;
  end
end
result = B{1};
marginal = pot_to_marginal(result);
[prob, mpe] = max(marginal.T);

% handle impossible cases
if ~(prob>0)
  mpe = [];    
  ll = -inf;
  %warning('evidence has zero probability')
  return
end

ll = log(prob);

% Do forward phase    
for ii=2:n
  marginal = pot_to_marginal(B{ii});
  mpeidx = [];
  for jj=order(1:length(mpe))
    assert(ismember(jj, marginal.domain)) %%% bug
    temp = find_equiv_posns(jj, marginal.domain);
    mpeidx = [mpeidx, temp] ;
    if isempty(temp)
      mpeidx = [mpeidx, Inf] ;
    end
  end
  [mpeidxsorted sortedtompe] = sort(mpeidx) ;
  
  % maximize the matrix obtained from assigning values from previous buckets.
  % this is done by building a string and using eval.
  
  kk=1;
  sargs = '(';
  for jj=1:length(marginal.domain)
    if (jj~=1)
      sargs = [sargs, ','];
    end
    if (mpeidxsorted(kk)==jj)
      sargs = [sargs, num2str(mpe(sortedtompe(kk)))];
      if (kk<length(mpe))
	kk = kk+1 ;
      end
    else
      sargs = [sargs, ':'];
    end
  end
  sargs = [sargs, ')'] ;   
  eval(['[val, loc] = max(marginal.T', sargs, ');'])        
  mpe = [mpe loc];
end     
[I,J] = sort(order);
mpe = mpe(J);



%%%%%%%%%

function b = bucket_num(domain, order)

b = max(find_equiv_posns(domain, order));

