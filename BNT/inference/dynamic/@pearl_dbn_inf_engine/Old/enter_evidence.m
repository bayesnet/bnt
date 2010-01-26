function [engine, loglik] = enter_evidence(engine, evidence, filter)
% ENTER_EVIDENCE Add the specified evidence to the network (pearl_dbn)
% [engine, loglik] = enter_evidence(engine, evidence, filter)
%
% evidence{i,t} = [] if if X(i,t) is hidden, and otherwise contains its observed value (scalar or column vector)
% If filter = 1, we do filtering, otherwise smoothing (default).

if nargin < 3, filter = 0; end

[ss T] = size(evidence);
bnet = bnet_from_engine(engine);
bnet2 = dbn_to_bnet(bnet, T);
ns = bnet2.node_sizes;
hnodes = mysetdiff(1:ss, engine.onodes);
hnodes = hnodes(:)';

[engine.parent_index, engine.child_index] = mk_pearl_msg_indices(bnet2);

msg = init_msgs(bnet2.dag, ns, evidence);
msg = init_ev_msgs(engine, evidence, msg);

niter = 1;
for iter=1:niter
  % FORWARD
  for t=1:T
    % update pi
    for i=1:ss %hnodes
      n = i + (t-1)*ss;
      ps = parents(bnet2.dag, n);
      if t==1
	e = bnet.equiv_class(i,1);
      else
	e = bnet.equiv_class(i,2);
      end
      msg{n}.pi = compute_pi(bnet.CPD{e}, n, ps, msg);
      %msg{n}.pi = normalise(msg{n}.pi(:) .* msg{n}.lambda_from_self(:));
    end
    % send pi msg to children
    for i=1:ss % hnodes
      n = i + (t-1)*ss;
      cs = children(bnet2.dag, n);
      for c=cs(:)'
	j = engine.parent_index{c}(n); % n is c's j'th parent
	msg{c}.pi_from_parent{j} = normalise(compute_pi_msg(n, cs, msg, c, ns));
      end
    end
  end

  if filter
    disp('skipping smoothing');
    break;
  end
    
  % BACKWARD
  for t=T:-1:1
    % update lambda
    for i=1:ss % hnodes
      n = i + (t-1)*ss;
      cs = children(bnet2.dag, n);
      msg{n}.lambda = compute_lambda(n, cs, msg, ns);
    end
    % send lambda msgs to parents
    for i=1:ss % hnodes
      n = i + (t-1)*ss;
      ps = parents(bnet2.dag, n);
      for p=ps(:)'
	j = engine.child_index{p}(n); % n is p's j'th child
	if t > 1
	  e = bnet.equiv_class(i, 2);
	else
	  e = bnet.equiv_class(i, 1);
	end
	msg{p}.lambda_from_child{j} = normalise(compute_lambda_msg(bnet.CPD{e}, n, ps, msg, p));
      end 
    end
  end
  
end


engine.marginal = cell(ss,T);
lik = zeros(1,ss*T);
for t=1:T
  for i=1:ss
    n = i + (t-1)*ss;
    [bel, lik(n)] = normalise(msg{n}.pi .* msg{n}.lambda);     
    engine.marginal{i,t} = bel;
  end
end

engine.evidence = evidence; % needed by marginal_nodes and marginal_family
engine.msg = msg;  % needed by marginal_family
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

lam = msg{n}.lambda_from_self(:);
%lam = ones(ns(n), 1);
for i=1:length(cs)
  c = cs(i);
  if c ~= except
    lam = lam .* msg{n}.lambda_from_child{i};
  end
end   

