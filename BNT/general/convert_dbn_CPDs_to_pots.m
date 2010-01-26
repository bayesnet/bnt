function CPDpot = convert_dbn_CPDs_to_pots(bnet, evidence, pot_type, softCPDpot)
% CONVERT_DBN_CPDS_TO_POTS Convert CPDs of (possibly instantiated) DBN nodes to potentials
% CPDpot = convert_dbn_CPDs_to_pots(bnet, evidence, pot_type, softCPDpot)
%
% CPDpot{n,t} is a potential containing P(n,t|pa(n,t), ev)
% softCPDpot{n,t} is a potential containing P(n,t|pa(n,t), ev) insted of using n's CPD

[ss T] = size(evidence);

if nargin < 4, softCPDpot = cell(ss,T); end
CPDpot = softCPDpot;

% Convert CPDs of instantiated nodes to potential form
t = 1;
for n=1:ss
  fam = family(bnet.dag, n);
  e = bnet.equiv_class(n, 1);
  if isempty(softCPDpot{n,t})
    CPDpot{n,t} = convert_to_pot(bnet.CPD{e}, pot_type, fam(:), evidence(:,1));
  end
end
for n=1:ss
  fam = family(bnet.dag, n, 2);
  e = bnet.equiv_class(n, 2);
  for t=2:T
    if isempty(softCPDpot{n,t})
      CPDpot{n,t} = convert_to_pot(bnet.CPD{e}, pot_type, fam(:), evidence(:,t-1:t));
    end
  end       
end
