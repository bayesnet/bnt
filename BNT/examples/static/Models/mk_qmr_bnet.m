function bnet = mk_qmr_bnet(G, inhibit, leak, prior, tabular_findings, onodes)
% MK_QMR_BNET Make a QMR model
% bnet = mk_qmr_bnet(G, inhibit, leak, prior)
%
% G(i,j) = 1 iff there is an arc from disease i to finding j
% inhibit(i,j) = inhibition probability on i->j arc
% leak(j) = inhibition prob. on leak->j arc
% prior(i) = prob. disease i is on
% tabular_findings = 1 means multinomial leaves (ignores leak/inhibit params)
%   = 0 means noisy-OR leaves (default = 0)

if nargin < 5, tabular_findings = 0; end

[Ndiseases Nfindings] = size(inhibit);
N = Ndiseases + Nfindings;
finding_node = Ndiseases+1:N;
ns = 2*ones(1,N);
dag = zeros(N,N);
dag(1:Ndiseases, finding_node) = G;
if nargin < 6, onodes = finding_node; end
bnet = mk_bnet(dag, ns, 'observed', onodes);

for d=1:Ndiseases
  CPT = [1-prior(d) prior(d)];
  bnet.CPD{d} = tabular_CPD(bnet, d, CPT');
end

for i=1:Nfindings
  fnode = finding_node(i);
  ps = parents(G, i);
  if tabular_findings
    bnet.CPD{fnode} = tabular_CPD(bnet, fnode); 
  else
    bnet.CPD{fnode} = noisyor_CPD(bnet, fnode, leak(i), inhibit(ps, i));
  end
end





