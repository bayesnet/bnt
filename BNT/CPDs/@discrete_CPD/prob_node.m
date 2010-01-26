function [P, p] = prob_node(CPD, self_ev, pev)
% PROB_NODE Compute prod_m P(x(i,m)| x(pi_i,m), theta_i) for node i (discrete)
% [P, p] = prob_node(CPD, self_ev, pev)
%
% self_ev(m) is the evidence on this node in case m.
% pev(i,m) is the evidence on the i'th parent in case m (if there are any parents).
% (These may also be cell arrays.)
%
% p(m) = P(x(i,m)| x(pi_i,m), theta_i) 
% P = prod p(m)

if iscell(self_ev), usecell = 1; else usecell = 0; end

ncases = length(self_ev);
sz = dom_sizes(CPD);

nparents = length(sz)-1;
if nparents == 0
  assert(isempty(pev));
else
  assert(isequal(size(pev), [nparents ncases]));
end

n = length(sz);
dom = 1:n;
p = zeros(1, ncases);
if isa(CPD, 'tabular_CPD')
  % speed up by looking up CPT using index Zhang Yimin  2001-12-31
  if usecell
    if nparents == 0
      data = [cell2num(self_ev)]; 
    else
      data = [cell2num(pev); cell2num(self_ev)]; 
    end
  else
    if nparents == 0
      data = [self_ev];
    else
      data = [pev; self_ev];
    end
  end
  
  indices = subv2ind(sz, data'); % each row of data' is a case 
  
  CPT=CPD_to_CPT(CPD);
  p = CPT(indices);
  
  %get the prob list
  %cpt_size = prod(sz);
  %prob_list=reshape(CPT, cpt_size, 1);
  %for m=1:ncases  %here we assume we get evidence for node and all its parents
  %  idx=indices(m);
  %  p(m)=prob_list(idx); 
  %end
  
else % eg. softmax
  
  for m=1:ncases
    if usecell
      if nparents == 0
	evidence = {self_ev{m}};
      else
	evidence = cell(1,n);
	evidence(1:n-1) = pev(:,m);
	evidence(n) = self_ev(m);
      end
    else
      if nparents == 0
	evidence = num2cell(self_ev(m));
      else
	evidence = num2cell([pev(:,m)', self_ev(m)]);
      end
    end
    T = convert_to_table(CPD, dom, evidence);
    p(m) = T;
  end
end
  
P = prod(p);


