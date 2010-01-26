function bnet = mk_asia_bnet(CPD_type, p, arity)
% MK_ASIA_BNET Make the 'Asia' bayes net.
%
% BNET = MK_ASIA_BNET uses the parameters specified on p21 of Cowell et al, 
% "Probabilistic networks and expert systems", Springer Verlag 1999.
% 
% BNET = MK_ASIA_BNET('cpt', p) uses random parameters drawn from a Dirichlet(p,p,...)
% distribution. If p << 1, this is nearly deterministic; if p >> 1, this is nearly uniform.
% 
% BNET = MK_ASIA_BNET('bool') makes each CPT a random boolean function.
%
% BNET = MK_ASIA_BNET('gauss') makes each CPT a random linear Gaussian distribution.
%
% BNET = MK_ASIA_BNET('orig') is the same as MK_ASIA_BNET.
%
% BNET = MK_ASIA_BNET('cpt', p, arity) can specify non-binary nodes.


if nargin == 0, CPD_type = 'orig'; end
if nargin < 3, arity = 2; end

Smoking = 1;
Bronchitis = 2;
LungCancer = 3;
VisitToAsia = 4;
TB = 5;
TBorCancer = 6;
Dys = 7;
Xray = 8;

n = 8;
dag = zeros(n);
dag(Smoking, [Bronchitis LungCancer]) = 1;
dag(Bronchitis, Dys) = 1;
dag(LungCancer, TBorCancer) = 1;
dag(VisitToAsia, TB) = 1;
dag(TB, TBorCancer) = 1;
dag(TBorCancer, [Dys Xray]) = 1;

ns = arity*ones(1,n);
if strcmp(CPD_type, 'gauss')
  dnodes = [];
else
  dnodes = 1:n;
end
bnet = mk_bnet(dag, ns, 'discrete', dnodes);

switch CPD_type
  case 'orig', 
    % true is 2, false is 1
    bnet.CPD{VisitToAsia} = tabular_CPD(bnet, VisitToAsia, [0.99   0.01]);
    bnet.CPD{Bronchitis} = tabular_CPD(bnet, Bronchitis, [0.7 0.4   0.3 0.6]);
    % minka: bug fix
    bnet.CPD{Dys} = tabular_CPD(bnet, Dys, [0.9 0.2 0.3 0.1   0.1 0.8 0.7 0.9]);
    bnet.CPD{TBorCancer} = tabular_CPD(bnet, TBorCancer, [1 0 0 0   0 1 1 1]);
    % minka: bug fix
    bnet.CPD{LungCancer} = tabular_CPD(bnet, LungCancer, [0.99 0.9  0.01 0.1]);
    bnet.CPD{Smoking} = tabular_CPD(bnet, Smoking, [0.5 0.5]);
    bnet.CPD{TB} = tabular_CPD(bnet, TB, [0.99 0.95  0.01 0.05]);
    bnet.CPD{Xray} = tabular_CPD(bnet, Xray, [0.95 0.02  0.05 0.98]);
 case 'bool',
  for i=1:n
    bnet.CPD{i} = boolean_CPD(bnet, i, 'rnd');
  end
 case 'gauss',
  for i=1:n
    bnet.CPD{i} = gaussian_CPD(bnet, i, 'cov', 1*eye(ns(i)));
  end
 case 'cpt',
  for i=1:n
    bnet.CPD{i} = tabular_CPD(bnet, i, p);
  end
end

  
  
