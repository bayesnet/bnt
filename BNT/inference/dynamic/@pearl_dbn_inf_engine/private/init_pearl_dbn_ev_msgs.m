function msg = init_pearl_dbn_ev_msgs(bnet, evidence, engine)

[ss T] = size(evidence);
pot_type = 'd';

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
      n = i + (t-1)*ss;
      lam_msg = normalise(temp.T);
      j = engine.child_index{n}(c+(t-1)*ss);
      assert(j==1);
      msg{n}.lambda_from_child{j} = lam_msg;
    end
  end
end
