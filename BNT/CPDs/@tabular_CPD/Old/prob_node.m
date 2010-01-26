function p = prob_node(CPD, self_ev, pev)
% PROB_NODE Compute P(y|pa(y), theta) (tabular)
% p = prob_node(CPD, self_ev, pev)
%
% self_ev{m} is the evidence on this node in case m
% pev{i,m} is the evidence on the i'th parent in case m
% If there is a single case, self_ev can be a scalar instead of a cell array

ncases = size(pev, 2);

%assert(~any(isemptycell(pev))); % slow
%assert(~any(isemptycell(self_ev))); % slow

CPT = CPD_to_CPT(CPD);  
sz = mysize(CPT);
nparents = length(sz)-1;
assert(nparents == size(pev, 1));

if ncases==1 
  x = cat(1, pev{:});
  if iscell(y)
    y = self_ev{1};
  else
    y = self_ev;
  end
  switch nparents
   case 0, p = CPT(y);
   case 1, p = CPT(x(1), y);
   case 2, p = CPT(x(1), x(2), y);
   case 3, p = CPT(x(1), x(2), x(3), y);
   otherwise,
    ind = subv2ind(CPD.sizes, [x y]);
    p = CPT(ind);
  end
else
  x = num2cell(pev)'; % each row is a case
  y = cat(1, self_ev{:})';
  ind = subv2ind(CPD.sizes, [x y]);
  p = CPT(ind);
end     
