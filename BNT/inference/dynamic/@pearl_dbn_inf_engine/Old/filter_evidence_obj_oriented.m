function [marginal, msg, loglik] = filter_evidence_old(engine, evidence)
% [marginal, msg, loglik] = filter_evidence(engine, evidence) (pearl_dbn)

[ss T] = size(evidence);
bnet = bnet_from_engine(engine);
bnet2 = dbn_to_bnet(bnet, T);
ns = bnet2.node_sizes;
hnodes = mysetdiff(1:ss, engine.onodes);
hnodes = hnodes(:)';

[engine.parent_index, engine.child_index] = mk_pearl_msg_indices(bnet2);

msg = init_msgs(bnet2.dag, ns, evidence);
msg = init_ev_msgs(engine, evidence, msg);

verbose = 1;
if verbose, fprintf('\nold filtering\n'); end

for t=1:T
  % update pi
  for i=hnodes
    n = i + (t-1)*ss;
    ps = parents(bnet2.dag, n);
    if t==1
      e = bnet.equiv_class(i,1);
    else
      e = bnet.equiv_class(i,2);
    end
    msg{n}.pi = compute_pi(bnet.CPD{e}, n, ps, msg);
    %if verbose, fprintf('%d computes pi\n', n); disp(msg{n}.pi); end
    msg{n}.pi = normalise(msg{n}.pi(:) .* msg{n}.lambda_from_self(:));
    if verbose, fprintf('%d recomputes pi\n', n); disp(msg{n}.pi); end
  end
  % send pi msg to children
  for i=hnodes
    n = i + (t-1)*ss;
    cs = children(bnet2.dag, n);
    for c=cs(:)'
      j = engine.parent_index{c}(n); % n is c's j'th parent
      pi_msg = normalise(compute_pi_msg(n, cs, msg, c, ns));
      msg{c}.pi_from_parent{j} = pi_msg;
      if verbose, fprintf('%d sends pi to %d\n', n,c); disp(pi_msg); end
    end
  end
end


marginal = cell(ss,T);
lik = zeros(1,ss*T);
for t=1:T
  for i=1:ss
    n = i + (t-1)*ss;
    %[bel, lik(n)] = normalise(msg{n}.pi .* msg{n}.lambda);     
    [bel, lik(n)] = normalise(msg{n}.pi);
    marginal{i,t} = bel;
  end
end

loglik = sum(log(lik));



%%%%%%%

function lambda = compute_lambda(n, cs, msg, ns)
% Pearl p183 eq 4.50
lambda = prod_lambda_msgs(n, cs, msg, ns);

%%%%%%%

function pi_msg = compute_pi_msg(n, cs, msg, c, ns)
% Pearl p183 eq 4.53 and 4.51
pi_msg = msg{n}.pi .* prod_lambda_msgs(n, cs, msg, ns, c);

%%%%%%%%%

function lam = prod_lambda_msgs(n, cs, msg, ns, except)

if nargin < 5, except = -1; end

%lam = msg{n}.lambda_from_self(:);
lam = ones(ns(n), 1);
for i=1:length(cs)
  c = cs(i);
  if c ~= except
    lam = lam .* msg{n}.lambda_from_child{i};
  end
end   


%%%%%%%%%%%

function msg = init_msgs(dag, ns, evidence)
% INIT_MSGS Initialize the lambda/pi message and state vectors (pearl_dbn)
% msg =  init_msgs(dag, ns, evidence)
%
% We assume all the hidden nodes are discrete.

N = length(dag);
msg = cell(1,N);
observed = ~isemptycell(evidence(:));

for n=1:N
  ps = parents(dag, n);
  msg{n}.pi_from_parent = cell(1, length(ps));
  for i=1:length(ps)
    p = ps(i);
    msg{n}.pi_from_parent{i} = ones(ns(p), 1);
  end
  
  cs = children(dag, n);
  msg{n}.lambda_from_child = cell(1, length(cs));
  for i=1:length(cs)
    c = cs(i);
    msg{n}.lambda_from_child{i} = ones(ns(n), 1);
  end

  msg{n}.lambda = ones(ns(n), 1);
  msg{n}.pi = ones(ns(n), 1);
  
  msg{n}.lambda_from_self = ones(ns(n), 1);
end


%%%%%%%%%

function msg = init_ev_msgs(engine, evidence, msg)
% Initialize the lambdas with any evidence

[ss T] = size(evidence);
bnet = bnet_from_engine(engine);
pot_type = 'd';
t = 1;
hnodes = mysetdiff(1:ss, engine.onodes);
for i=hnodes(:)'
  c = engine.obschild(i);
  if c > 0
    fam = family(bnet.dag, c);
    e = bnet.equiv_class(c, 1);
    CPDpot = CPD_to_pot(pot_type, bnet.CPD{e}, fam, bnet.node_sizes(:), bnet.cnodes(:), evidence(:,1));
    temp = pot_to_marginal(CPDpot);
    n = i;
    msg{n}.lambda_from_self = temp.T;
  end
end
for t=2:T
  for i=hnodes(:)'
    c = engine.obschild(i);
    if c > 0 
      fam = family(bnet.dag, c, 2);
      e = bnet.equiv_class(c, 2);
      CPDpot = CPD_to_pot(pot_type, bnet.CPD{e}, fam, bnet.node_sizes(:), bnet.cnodes(:), evidence(:,t-1:t));
      temp = pot_to_marginal(CPDpot);
      n = i + (t-1)*ss;
      msg{n}.lambda_from_self = temp.T;
    end
  end
end       
