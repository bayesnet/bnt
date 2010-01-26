function bnet = mk_ideker_bnet(CPD_type, p)
% MK_IDEKER_BNET Make the Bayes net in the PSB'00 paper by Ideker, Thorsson and Karp.
%
% BNET = MK_IDEKER_BNET uses the boolean functions specified in the paper 
% "Discovery of regulatory interactions through perturbation: inference and experimental design",
% Pacific Symp. on Biocomputing, 2000.
% 
% BNET = MK_IDEKER_BNET('root') uses the above boolean functions, but puts a uniform
% distribution on the root nodes.
%
% BNET = MK_IDEKER_BNET('cpt', p) uses random parameters drawn from a Dirichlet(p,p,...)
% distribution. If p << 1, this is nearly deterministic; if p >> 1, this is nearly uniform.
% 
% BNET = MK_IDEKER_BNET('bool') makes each CPT a random boolean function.
%
% BNET = MK_IDEKER_BNET('orig') is the same as MK_IDEKER_BNET.


if nargin == 0
  CPD_type = 'orig';
end

n = 4;
dag = zeros(n);
dag(1,3)=1;
dag(2,[3 4])=1;
dag(3,4)=1;
ns = 2*ones(1,n);
bnet = mk_bnet(dag, ns);

switch CPD_type
 case 'orig',
  bnet.CPD{1} = tabular_CPD(bnet, 1, [0 1]);
  bnet.CPD{2} = tabular_CPD(bnet, 2, [0 1]);
  bnet.CPD{3} = boolean_CPD(bnet, 3, 'inline', inline('x(1) & x(2)'));
  bnet.CPD{4} = boolean_CPD(bnet, 4, 'inline', inline('x(1) & ~x(2)'));
 case 'root',
  bnet.CPD{1} = tabular_CPD(bnet, 1, [0.5 0.5]);
  bnet.CPD{2} = tabular_CPD(bnet, 2, [0.5 0.5]);
  bnet.CPD{3} = boolean_CPD(bnet, 3, 'inline', inline('x(1) & x(2)'));
  bnet.CPD{4} = boolean_CPD(bnet, 4, 'inline', inline('x(1) & ~x(2)'));
 case 'bool',
  for i=1:n
    bnet.CPD{i} = boolean_CPD(bnet, i, 'rnd');
  end
 case 'cpt',
  for i=1:n
    bnet.CPD{i} = tabular_CPD(bnet, i, p);
  end
 otherwise,
  error(['unknown type ' CPD_type]);
end
