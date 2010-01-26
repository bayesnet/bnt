function [marginal, msg, loglik] = filter_evidence(engine, evidence)

error('broken');

[ss T] = size(evidence);
bnet = bnet_from_engine(engine);
onodes = engine.onodes;
hnodes = mysetdiff(1:ss, onodes);
hnodes = hnodes(:)';

ns = bnet.node_sizes(:);
onodes2 = [onodes(:); onodes(:)+ss];
ns(onodes2) = 1;
	   
verbose = 1;
if verbose, fprintf('\nnew filtering\n'); end
  
pot_type = 'd';
niter = engine.max_iter;

% msg(i1,t1,i2,j2) (i1,t1) -> (i2,t2)
%lambda_msg = cell(ss,T,ss,T);
%pi_msg = cell(ss,T,ss,T);

% intra_lambda_msg(i,j,t) (i,t) -> (j,t), i is child
% inter_lambda_msg(i,j,t) (i,t+1) -> (j,t), i is child
% inter_pi_msg(i,j,t) (i,t-1) -> (j,t), i is parent
intra_lambda_msg = cell(ss,ss,T);
inter_lambda_msg = cell(ss,ss,T);
inter_pi_msg = cell(ss,ss,T);

lambda = cell(ss,T);
pi = cell(ss,T);

for t=1:T
  for i=1:ss
    lambda{i,t} = ones(ns(i), 1);
    pi{i,t} = ones(ns(i), 1);
    
    cs = children(bnet.intra, i);
    for c=cs(:)'
      intra_lambda_msg{c,i,t} = ones(ns(i),1);
    end
    
    cs = children(bnet.inter, i);
    for c=cs(:)'
      inter_lambda_msg{c,i,t} = ones(ns(i),1);
    end
    
    ps = parents(bnet.inter, i);
    for p=ps(:)'
      inter_pi_msg{p,i,t} = ones(ns(i), 1); % not used for t==1
    end
  end
end

% each hidden node absorbs lambda from its observed child (if any)
for t=1:T
  for i=hnodes
    c = engine.obschild(i);
    if c > 0
      if t==1
	fam = family(bnet.dag, c);
	e = bnet.equiv_class(c, 1);
	CPDpot = CPD_to_pot(pot_type, bnet.CPD{e}, fam, bnet.node_sizes(:), bnet.cnodes(:), evidence(:,1));
      else
	fam = family(bnet.dag, c, 2); % within 2 slice network
	e = bnet.equiv_class(c, 2);
	CPDpot = CPD_to_pot(pot_type, bnet.CPD{e}, fam, bnet.node_sizes(:), bnet.cnodes(:), evidence(:,t-1:t));
      end
      temp = pot_to_marginal(CPDpot);
      lam_msg = normalise(temp.T);
      %if verbose, fprintf('(%d,%d) sends lambda to (%d,%d)\n', c,t, i,t); disp(lam_msg); end
      intra_lambda_msg{c,i,t} = lam_msg;
    end
  end
end

% FORWARD
for t=1:T
  % update pi
  for i=hnodes
    if t==1
      e = bnet.equiv_class(i,1);
      temp = struct(bnet.CPD{e});
      pi{i,t} = temp.CPT;
    else
      e = bnet.equiv_class(i,2);
      temp = struct(bnet.CPD{e});
      ps = parents(bnet.inter, i);
      dom = [ps i+ss];
      pot = dpot(dom, ns(dom), temp.CPT);
      for p=ps(:)'
	temp = dpot(p, ns(p), inter_pi_msg{p,i,t});
	pot = multiply_by_pot(pot, temp);
      end
      pot = marginalize_pot(pot, i+ss);
      temp = pot_to_marginal(pot);
      pi{i,t} = temp.T;
      %if verbose, fprintf('(%d,%d) computes pi\n', i,t); disp(pi{i,t}); end
    end
    
    c = engine.obschild(i);
    if c > 0
      pi{i,t} = normalise(pi{i,t} .* intra_lambda_msg{c,i,t});
    end
    %if verbose, fprintf('(%d,%d) recomputes pi\n', i,t); disp(pi{i,t}); end
    if verbose, fprintf('%d recomputes pi\n', i+(t-1)*ss); disp(pi{i,t}); end
  end
  
  % send pi msg to children 
  for i=hnodes
    cs = children(bnet.inter, i);
    for c=cs(:)'
      pot = pi{i,t};
      for k=cs(:)'
	if k ~= c
	  pot = pot .* inter_lambda_msg{k,i,t};
	end
      end
      cs2 = children(bnet.intra, i);
      for k=cs2(:)'
	pot = pot .* intra_lambda_msg{k,i,t};
      end
      pot = normalise(pot);
      %if verbose, fprintf('(%d,%d) sends pi to (%d,%d)\n', i,t, c,t+1); disp(pot); end
      if verbose, fprintf('%d sends pi to %d\n', i+(t-1)*ss, c+t*ss); disp(pot); end
      inter_pi_msg{i,c,t+1} = pot;
    end
  end
end


marginal = cell(ss,T);
for t=1:T
  for i=hnodes
    %marginal{i,t} = normalise(pi{i,t} .* lambda{i,t});     
    marginal{i,t} = normalise(pi{i,t});
  end
end

loglik = 0;

msg.inter_pi_msg = inter_pi_msg;
msg.inter_lambda_msg = inter_lambda_msg;
msg.intra_lambda_msg = intra_lambda_msg;
