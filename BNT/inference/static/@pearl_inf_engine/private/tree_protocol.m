function msg = tree_protocol(engine, evidence, msg)

bnet = bnet_from_engine(engine);
N = length(bnet.dag);

% Send messages from leaves to root
for i=1:N-1
  n = engine.postorder(i);
  above = parents(engine.adj_mat, n);
  msg = send_msgs_to_some_neighbors(n, msg, above, bnet, engine.child_index, engine.parent_index, ...
				    engine.msg_type, evidence);
end

% Process root
n = engine.root;
cs = children(bnet.dag, n);
%msg{n}.lambda = compute_lambda(n, cs, msg, engine.msg_type);
msg{n}.lambda = prod_lambda_msgs(n, cs, msg, engine.msg_type);
ps = parents(bnet.dag, n);
msg{n}.pi = CPD_to_pi(bnet.CPD{bnet.equiv_class(n)}, engine.msg_type, n, ps, msg, evidence);

% Send messages from root to leaves
for i=1:N
  n = engine.preorder(i);
  below = children(engine.adj_mat, n);
  msg = send_msgs_to_some_neighbors(n, msg, below, bnet, engine.child_index, engine.parent_index, ...
				    engine.msg_type, evidence);
end

  
%%%%%%%%%%

function msg = send_msgs_to_some_neighbors(n, msg, valid_nbrs, bnet, child_index, parent_index, ...
					   msg_type, evidence)

verbose = 0;

ns = bnet.node_sizes;
dag = bnet.dag;
e = bnet.equiv_class(n);
CPD = bnet.CPD{e};


cs = children(dag, n);
%msg{n}.lambda = compute_lambda(n, cs, msg);
msg{n}.lambda = prod_lambda_msgs(n, cs, msg, msg_type);
if verbose, fprintf('%d computes lambda\n', n); display(msg{n}.lambda); end

ps = parents(dag, n);
msg{n}.pi = CPD_to_pi(CPD, msg_type, n, ps, msg, evidence);
if verbose, fprintf('%d computes pi\n', n); display(msg{n}.pi); end

ps2 = myintersect(parents(dag, n), valid_nbrs);
for p=ps2(:)'
  lam_msg = CPD_to_lambda_msg(CPD, msg_type, n, ps, msg, p, evidence);
  j = child_index{p}(n); % n is p's j'th child
  msg{p}.lambda_from_child{j} = lam_msg;
  if verbose, fprintf('%d sends lambda to %d\n', n, p); display(lam_msg); end
end

cs2 = myintersect(cs, valid_nbrs);
for c=cs2(:)'
  %pi_msg = compute_pi_msg(n, cs, msg, c);
  pi_msg = compute_bel(msg_type, msg{n}.pi, prod_lambda_msgs(n, cs, msg, msg_type, c));
  j = parent_index{c}(n); % n is c's j'th parent
  msg{c}.pi_from_parent{j} = pi_msg;
  if verbose, fprintf('%d sends pi to %d\n', n, c); display(pi_msg); end
end


 
