% Make the following network (from Jensen (1996) p84 fig 4.17)
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

ns = [5 4 3 2 2 1 2 2 2]; % vector-valued nodes
%ns = ones(1,9); % scalar nodes
dnodes = [];

bnet = mk_bnet(dag, ns, 'discrete', []);
rand('state', 0);
randn('state', 0);
for i=1:N
  bnet.CPD{i} = gaussian_CPD(bnet, i);
end

clear engine;
engine{1} = gaussian_inf_engine(bnet);
engine{2} = jtree_inf_engine(bnet);

[err, time] = cmp_inference_static(bnet, engine);

