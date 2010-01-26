function [msg, niter] = parallel_protocol(engine, evidence, msg)

bnet = bnet_from_engine(engine);
N = length(bnet.dag);
ns = bnet.node_sizes(:);

if ~isempty(engine.filename)
  fid = fopen(engine.filename, 'w');
  if fid == 0
    error(['could not open ' engine.filename ' for writing'])
  end
else
  fid = [];
end

converged = 0;
iter = 1;
hidden = find(isemptycell(evidence));
bel = cell(1,N);
old_bel = cell(1,N);
%nodes = mysetdiff(1:N, engine.disconnected_nodes);
nodes = find(~engine.disconnected_nodes_bitv);
while ~converged & (iter <= engine.max_iter)
  % Everybody updates their state in parallel
  for n=nodes(:)'
    cs_msg = children(engine.msg_dag, n);
    %msg{n}.lambda = compute_lambda(n, cs, msg);
    msg{n}.lambda = prod_lambda_msgs(n, cs_msg, msg, engine.msg_type);
    ps_orig = parents(bnet.dag, n);
    msg{n}.pi = CPD_to_pi(bnet.CPD{bnet.equiv_class(n)}, engine.msg_type, n, ps_orig, msg, evidence);
  end
  
  changed = 0;
  if ~isempty(fid)
    fprintf(fid, 'ITERATION %d\n', iter);
  end
  for n=hidden(:)' % this will not contain any disconnected nodes
    old_bel{n} = bel{n};
    bel{n}  = compute_bel(engine.msg_type, msg{n}.pi, msg{n}.lambda);
    if ~isempty(fid)
      fprintf(fid, 'node %d: %s\n', n, bel_to_str(bel{n}, engine.msg_type));
    end
    if engine.storebel
      engine.bel{n,iter} = bel{n};
    end
    if (iter == 1) | ~approxeq_bel(bel{n}, old_bel{n}, engine.tol, engine.msg_type)
      changed = 1;
    end
  end
  %converged = ~changed;
  converged = ~changed & (iter > 1);  % Sonia Leach changed this

  if ~converged
    % Everybody sends to all their neighbors in parallel
    for n=nodes(:)'
      % lambda msgs to parents
      ps_msg = parents(engine.msg_dag, n);
      ps_orig = parents(bnet.dag, n);
      for p=ps_msg(:)'
	j = engine.child_index{p}(n); % n is p's j'th child
	old_msg = msg{p}.lambda_from_child{j}(:);
	new_msg = CPD_to_lambda_msg(bnet.CPD{bnet.equiv_class(n)}, engine.msg_type, n, ps_orig, ...
				    msg, p, evidence);
	lam_msg = convex_combination_msg(old_msg, new_msg, engine.momentum, engine.msg_type);
	msg{p}.lambda_from_child{j} = lam_msg;
      end 

      % pi msgs to children
      cs_msg = children(engine.msg_dag, n);
      for c=cs_msg(:)'
	j = engine.parent_index{c}(n); % n is c's j'th parent
	old_msg = msg{c}.pi_from_parent{j}(:);
	%new_msg = compute_pi_msg(n, cs, msg, c));
	new_msg = compute_bel(engine.msg_type, msg{n}.pi, prod_lambda_msgs(n, cs_msg, msg, engine.msg_type, c));
	pi_msg = convex_combination_msg(old_msg, new_msg, engine.momentum, engine.msg_type);
	msg{c}.pi_from_parent{j} = pi_msg;
      end
    end
    iter = iter + 1;
  end
end

if fid > 0, fclose(fid); end
%niter = iter - 1;
niter = iter;

%%%%%%%%%%

function str = bel_to_str(bel, type)

switch type
 case 'd', str = sprintf('%9.4f ', bel(:)');
 case 'g', str = sprintf('%9.4f ', bel.mu(:)');
end


%%%%%%%

function a = approxeq_bel(bel1, bel2, tol, type)

switch type
 case 'd', a = approxeq(bel1, bel2, tol);
 case 'g', a = approxeq(bel1.mu, bel2.mu, tol) & approxeq(bel1.Sigma, bel2.Sigma, tol);
end


%%%%%%%

function msg = convex_combination_msg(old_msg, new_msg, old_weight, type)

switch type
 case 'd', msg = old_weight * old_msg + (1-old_weight)*new_msg;
 case 'g', msg = new_msg;
end
