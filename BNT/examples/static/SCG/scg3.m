% Compare various inference engines on the following network (from Jensen (1996) p84 fig 4.17)
%    1
%  / | \
% 2  3  4
% |  |  |
% 5  6  7
%  \/ \/
%  8   9
% where all arcs point downwards

N = 9;
dag = zeros(N,N);
dag(1,2)=1; dag(1,3)=1; dag(1,4)=1;
dag(2,5)=1; dag(3,6)=1; dag(4,7)=1;
dag(5,8)=1; dag(6,8)=1; dag(6,9)=1; dag(7,9) = 1;

gauss = 1;
if gauss
  ns = ones(1,N); % scalar nodes
  ns(1) = 2;
  ns(9) = 3;
  dnodes = [];
else
  ns = 2*ones(1,N); % binary nodes
  dnodes = 1:N;
end

bnet = mk_bnet(dag, ns, 'discrete', dnodes);
% use random params
for i=1:N
  if gauss
    bnet.CPD{i} = gaussian_CPD(bnet, i);
  else
    bnet.CPD{i} = tabular_CPD(bnet, i);
  end
end

engines = {};
engines{1} = jtree_inf_engine(bnet);
engines{2} = stab_cond_gauss_inf_engine(bnet);

[err, time] = cmp_inference_static(bnet, engines);
