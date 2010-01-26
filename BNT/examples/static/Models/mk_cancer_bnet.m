function bnet = mk_cancer_bnet(CPD_type, p)
% MK_CANCER_BNET Make the 'Cancer' Bayes net.
%
% BNET = MK_CANCER_BNET uses the noisy-or parameters specified in Fig 4a of the UAI98 paper by
% Friedman, Murphy and Russell, "Learning the Structure of DPNs", p145.
%
% BNET = MK_CANCER_BNET('noisyor', p) makes each CPD a noisy-or, with probability p of
% suppression for each parent; leaks are turned off.
%
% BNET = MK_CANCER_BNET('cpt', p) uses random CPT parameters drawn from a Dirichlet(p,p,...)
% distribution. If p << 1, this is near deterministic; if p >> 1, this is near 1/k.
% p defaults to 1.0 (uniform distribution).
%
% BNET = MK_CANCER_BNET('bool') makes each CPT a random boolean function.
%
% In all cases, the root is set to a uniform distribution.

if nargin == 0
  rnd = 0;
else
  rnd = 1;
end

n = 5;
dag = zeros(n);
dag(1,[2 3]) = 1;
dag(2,4) = 1;
dag(3,4) = 1;
dag(4,5) = 1;

ns = 2*ones(1,n);
bnet = mk_bnet(dag, ns);
    
if ~rnd
  bnet.CPD{1} = tabular_CPD(bnet, 1, [0.5 0.5]);
  bnet.CPD{2} = noisyor_CPD(bnet, 2, 1.0, 1-0.9);
  bnet.CPD{3} = noisyor_CPD(bnet, 3, 1.0, 1-0.2);
  bnet.CPD{4} = noisyor_CPD(bnet, 4, 1.0, 1-[0.7 0.6]);
  bnet.CPD{5} = noisyor_CPD(bnet, 5, 1.0, 1-0.5);
else
  switch CPD_type
   case 'noisyor',
    for i=1:n
      ps = parents(dag, i);
      bnet.CPD{i} = noisyor_CPD(bnet, i, 1.0, p*ones(1,length(ps)));
    end
   case 'bool',
    for i=1:n
      bnet.CPD{i} = boolean_CPD(bnet, i, 'rnd');
    end
   case 'cpt',
    for i=1:n
      bnet.CPD{i} = tabular_CPD(bnet, i, p);
    end
   otherwise
    error(['bad CPD type ' CPD_type]);
  end
end
  
  
  
