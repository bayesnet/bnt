function [marginal, msg, loglik] = smooth_evidence(engine, evidence)
% [marginal, msg, loglik] = smooth_evidence(engine, evidence) (pearl_dbn)

disp('warning: pearl_dbn smoothing is broken');

[ss T] = size(evidence);
bnet = bnet_from_engine(engine);
bnet2 = dbn_to_bnet(bnet, T);
ns = bnet2.node_sizes;
hnodes = mysetdiff(1:ss, engine.onodes);
hnodes = hnodes(:)';

onodes2 = unroll_set(engine.onodes(:), ss, T);
onodes2 = onodes2(:)';

hnodes2 = unroll_set(hnodes(:), ss, T);
hnodes2 = hnodes2(:)';

[engine.parent_index, engine.child_index] = mk_pearl_msg_indices(bnet2);

msg = init_msgs(bnet2.dag, ns, evidence, bnet2.equiv_class, bnet2.CPD);

verbose = 0;
pot_type = 'd';
niter = 1;
for iter=1:niter
  % FORWARD
  for t=1:T
    if verbose, fprintf('t=%d\n', t); end

    % each hidden node absorbs lambda from its observed child (if any)
    for i=hnodes
      c = engine.obschild(i);
      if c > 0
	if t==1
	  fam = family(bnet.dag, c);
	  e = bnet.equiv_class(c, 1);
	  CPDpot = CPD_to_pot(pot_type, bnet.CPD{e}, fam, bnet.node_sizes(:), bnet.cnodes(:), evidence(:,1));
	else
	  fam = family(bnet.dag, 2); % within 2 slice network
	  e = bnet.equiv_class(c, 2);
	  CPDpot = CPD_to_pot(pot_type, bnet.CPD{e}, fam, bnet.node_sizes(:), bnet.cnodes(:), evidence(:,t-1:t));
	end
	temp = pot_to_marginal(CPDpot);
	n = i + (t-1)*ss;
	lam_msg = normalise(temp.T);
	j = engine.child_index{n}(c+(t-1)*ss);
	assert(j==1);
	msg{n}.lambda_from_child{j} = lam_msg;
	if verbose, fprintf('%d sends lambda to %d\n', c + (t-1)*ss, n); disp(lam_msg); end
      end
    end
    
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
      if verbose, fprintf('%d computes pi\n', n); disp(msg{n}.pi); end
    end
    
    % send pi msg to children in next slice
    for i=hnodes
      n = i + (t-1)*ss;
      %cs = myintersect(children(bnet2.dag, n), hnodes2);
      cs = children(bnet2.dag, n);
      for c=cs(:)'
	j = engine.parent_index{c}(n); % n is c's j'th parent
	pi_msg = normalise(compute_pi_msg(n, cs, msg, c, ns));
	msg{c}.pi_from_parent{j} = pi_msg;
	if verbose, fprintf('%d sends pi to %d\n', n, c); disp(pi_msg); end
      end
    end
  end

  % BACKWARD
  for t=T:-1:1
    if verbose, fprintf('t = %d\n', t); end

    % update lambda
    for i=hnodes
      n = i + (t-1)*ss;
      cs = children(bnet2.dag, n);
      msg{n}.lambda = compute_lambda(n, cs, msg, ns);
      if verbose, fprintf('%d computes lambda\n', n); disp(msg{n}.lambda); end
    end
    
    % send lambda msgs to hidden parents in prev slcie
    for i=hnodes
      n = i + (t-1)*ss;
      %ps = myintersect(parents(bnet2.dag, n), hnodes2);
      ps = parents(bnet2.dag, n);
      for p=ps(:)'
	j = engine.child_index{p}(n); % n is p's j'th child
	if t > 1
	  e = bnet.equiv_class(i, 2);
	else
	  e = bnet.equiv_class(i, 1);
	end
	lam_msg = normalise(compute_lambda_msg(bnet.CPD{e}, n, ps, msg, p));
	msg{p}.lambda_from_child{j} = lam_msg;
	if verbose, fprintf('%d sends lambda to %d\n', n, p); disp(lam_msg); end
      end 
    end
        
    % send pi msg to observed children 
    if 0
    for i=hnodes
      n = i + (t-1)*ss;
      cs = myintersect(children(bnet2.dag, n), onodes2);
      %cs = children(bnet2.dag, n);
      for c=cs(:)'
	j = engine.parent_index{c}(n); % n is c's j'th parent
	pi_msg = normalise(compute_pi_msg(n, cs, msg, c, ns));
	msg{c}.pi_from_parent{j} = pi_msg;
	if verbose, fprintf('%d sends pi to %d\n', n, c); disp(pi_msg); end
      end
    end
    end
    
  end
end


marginal = cell(ss,T);
lik = zeros(1,ss*T);
for t=1:T
  for i=hnodes
    n = i + (t-1)*ss;
    [bel, lik(n)] = normalise(msg{n}.pi .* msg{n}.lambda);     
    marginal{i,t} = bel;
  end
end

loglik = 0;
%loglik = sum(log(lik));



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


%%%%%%%%%

function msg = init_msgs(dag, ns, evidence, eclass, CPD)
% INIT_MSGS Initialize the lambda/pi message and state vectors (pearl_dbn)
% msg =  init_msgs(dag, ns, evidence)

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

  % Initialize the lambdas with any evidence
  if observed(n)
    v = evidence{n};
    msg{n}.lambda = zeros(ns(n), 1);
    msg{n}.lambda(v) = 1; % delta function
    msg{n}.lambda = [];
  end      
  
end

