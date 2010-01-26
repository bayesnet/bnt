function scg_unstable()

% the objective of this script is to test if the stable conditonal gaussian
% inference can handle the numerical instability problem described on
% page.151 of 'Probabilistic networks and expert system' by Cowell, Dawid, Lauritzen and
% Spiegelhalter, 1999.

A = 1; Y = 2;
n = 2;

ns = ones(1, n);
dnodes = [A];
cnodes = Y;
ns = [2 1];

dag = zeros(n);
dag(A, Y) = 1;

bnet = mk_bnet(dag, ns, dnodes);

bnet.CPD{A} = tabular_CPD(bnet, A, [0.5 0.5]'); 
bnet.CPD{Y} = gaussian_CPD(bnet, Y, 'mean', [0 1], 'cov', [1e-5 1e-6]);

evidence = cell(1, n);

pot_type = 'cg';
potYgivenA = convert_to_pot(bnet.CPD{Y}, pot_type, [A Y], evidence);
potA = convert_to_pot(bnet.CPD{A}, pot_type, A, evidence);
potYandA = multiply_by_pot(potYgivenA, potA);
potA2 = marginalize_pot(potYandA, A);

thresh = 1; % 0dp

[g,h,K] = extract_can(potA);
assert(approxeq(g(:)', [-0.693147 -0.693147], thresh))


[g,h,K] = extract_can(potYgivenA);
assert(approxeq(g(:)', [4.83752 -499994], thresh))
assert(approxeq(h(:)', [0 1e6]))
assert(approxeq(K(:)', [1e5 1e6]))

[g,h,K] = extract_can(potYandA);
assert(approxeq(g(:)', [4.14437 -499995], thresh))
assert(approxeq(h(:)', [0 1e6]))
assert(approxeq(K(:)', [1e5 1e6]))


[g,h,K] = extract_can(potA2);
%assert(approxeq(g(:)', [-0.69315 -1]))
g
assert(approxeq(g(:)', [-0.69315 -0.69315]))



if 0
pot_type = 'scg';
spotYgivenA = convert_to_pot(bnet.CPD{Y}, pot_type, [A Y], evidence);
spotA = convert_to_pot(bnet.CPD{A}, pot_type, A, evidence);
spotYandA = direct_combine_pots(spotYgivenA, spotA); 
spotA2 = marginalize_pot(spotYandA, A);

spotA=struct(spotA);
spotA2=struct(spotA2);
for i=1:2
  assert(approxeq(spotA2.scgpotc{i}.p, spotA.scgpotc{i}.p))
  assert(approxeq(spotA2.scgpotc{i}.A, spotA.scgpotc{i}.A))
  assert(approxeq(spotA2.scgpotc{i}.B, spotA.scgpotc{i}.B))
  assert(approxeq(spotA2.scgpotc{i}.C, spotA.scgpotc{i}.C))
end

end


%%%%%%%%%%%

function [g,h,K] = extract_can(pot)

pot = struct(pot);
D = length(pot.can);
g = zeros(1, D);
h = zeros(1, D);
K = zeros(1, D);
for i=1:D
  S = struct(pot.can{i});
  g(i) = S.g;
  if length(S.h) > 0
    h(i) = S.h;
    K(i) = S.K;
  end
end
