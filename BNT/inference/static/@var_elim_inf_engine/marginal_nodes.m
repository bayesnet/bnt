function [marginal, loglik] = marginal_nodes(engine, query, add_ev)
% MARGINAL_NODES Compute the marginal on the specified query nodes (var_elim)
% [marginal, loglik] = marginal_nodes(engine, query)

if nargin < 3, add_ev = 0; end

assert(length(query)>=1);

evidence = engine.evidence;

bnet = bnet_from_engine(engine);
ns = bnet.node_sizes;
n = length(bnet.dag);

onodes = find(~isemptycell(evidence));
hnodes = find(isemptycell(evidence));
pot_type = determine_pot_type(bnet, onodes);

% Fold the evidence into the CPTs - this could be done in 'enter_evidence'
CPT = cell(1,n);
for i=1:n
  fam = family(bnet.dag, i);
  CPT{i} = convert_to_pot(bnet.CPD{bnet.equiv_class(i)}, pot_type, fam(:), evidence);
end



sum_over = mysetdiff(1:n, query);
order = [query sum_over]; % no attempt to optimize this

% Initialize the buckets with the product of the CPTs assigned to them
B = cell(1,n+1); 
for b=1:n+1
  B{b} = mk_initial_pot(pot_type, [], [], [], []);
end
for i=1:n
  b = bucket_num(domain_pot(CPT{i}), order);
  B{b} = multiply_pots(B{b}, CPT{i});
end

% Do the marginalization
sum_over = sum_over(length(sum_over):-1:1); % reverse
for i=sum_over(:)'
  % summing over variable i which occurs in bucket j
  j = bucket_num(i, order);
  rest = mysetdiff(domain_pot(B{j}), i);
  % minka
  if ~isempty(rest)
    temp = marginalize_pot(B{j}, rest);
    b = bucket_num(domain_pot(temp), order);
    %fprintf('summing over bucket %d (var %d), putting result into bucket %d\n', j, i, b);
    B{b} = multiply_pots(B{b}, temp);
  end
end

% Combine all the remaining buckets into one
result = B{1};
for i=2:length(query)
  if ~isempty(domain_pot(B{i}))
    result = multiply_pots(result, B{i});
  end
end
[result, loglik] = normalize_pot(result);


marginal = pot_to_marginal(result);
% minka: from jtree_inf_engine
if add_ev
  bnet = bnet_from_engine(engine);
  %marginal = add_ev_to_dmarginal(marginal, engine.evidence, bnet.node_sizes);
  marginal = add_evidence_to_gmarginal(marginal, engine.evidence, bnet.node_sizes, bnet.cnodes);
end

%%%%%%%%%

function b = bucket_num(domain, order)

b = max(find_equiv_posns(domain, order));

