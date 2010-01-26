% Find out how big the cliques are in an HHMM as a function of depth
% (This is how we get the complexity bound of O(D K^{1.5D}).)

if 0
Qsize = [];
Fsize = [];
Nclqs = [];
end

ds = 1:15;

for d = ds
  allQ = 1;
  [intra, inter, Qnodes, Fnodes, Onode] = mk_hhmm_topo(d, allQ);
  
  N = length(intra);
  ns = 2*ones(1,N);
  
  bnet = mk_dbn(intra, inter, ns);
  for i=1:N
    bnet.CPD{i} = tabular_CPD(bnet, i);
  end
  
  if 0
    T = 5;
    dag = unroll_dbn_topology(intra, inter, T);
    engine = jtree_unrolled_dbn_inf_engine(bnet, T, 'constrained', 1);
    S = struct(engine);
    S1 = struct(S.sub_engine);
  end
  
  engine = jtree_dbn_inf_engine(bnet);
  S = struct(engine);
  J = S.jtree_struct;
  
  ss = 2*d+1;
  Qnodes2 = Qnodes + ss;
  QQnodes = [Qnodes Qnodes2];
  
  % find out how many Q nodes in each clique, and how many F nodes
  C = length(J.cliques);
  Nclqs(d) = 0;
  for c=1:C
    Qsize(c,d) = length(myintersect(J.cliques{c}, QQnodes));
    Fsize(c,d) = length(myintersect(J.cliques{c}, Fnodes));
    if length(J.cliques{c}) > 1 % exclude observed leaves
      Nclqs(d) = Nclqs(d) + 1;
    end
  end
  %pred_max_Qsize(d) = ceil(d+(d+1)/2);
  pred_max_Qsize(d) = ceil(1.5*d);
  
  fprintf('d=%d\n', d);
  %fprintf('D=%d, max F = %d. max Q = %d, pred max Q = %d\n', ...
	%  D, max(Fsize), max(Qsize), ceil(D+(D+1)/2));
	     
  %histc(Qsize,1:max(Qsize)) % how many of each size?
end % next d


Q = 2;
pred_mass = ds.*(Q.^ds) + Q.^(ceil(1.5 * ds))
pred_mass2 = Q.^(ceil(1.5 * ds))

for d=ds
  mass(d) = 0;
  for c=1:C
    mass(d) = mass(d) + Q^Qsize(c,d);
  end
end
    

if 0
%plot(ds, max(Qsize), 'o-',  ds, pred_max_Qsize, '*--');
%plot(ds, max(Qsize), 'o-',  ds, 1.5*ds, '*--');
%plot(ds, mass, 'o-',  ds, pred_mass, '*--');
D = 15;
%plot(ds(1:D), mass(1:D), 'bo-',  ds(1:D), pred_mass(1:D), 'g*--', ds(1:D), pred_mass2(1:D), 'k+-.');
plot(ds(1:D), log(mass(1:D)), 'bo-',  ds(1:D), log(pred_mass(1:D)), 'g*--', ds(1:D), log(pred_mass2(1:D)), 'k+-.');

grid on
xlabel('depth of hierarchy')
title('max num Q nodes in any clique vs. depth')
legend('actual', 'predicted')

%previewfig(gcf, 'width', 3, 'height', 1.5, 'color', 'bw');
%exportfig(gcf, '/home/cs/murphyk/WP/ConferencePapers/HHMM/clqsize2.eps', ...
%          'width', 3, 'height', 1.5, 'color', 'bw');   

end


if 0
for d=ds
  effnumclqs(d) = length(find(Qsize(:,d)>0));
end
ds = 1:10;
Qs = 2:10;
maxC = size(Qsize, 1);
cost = [];
cost_bound = [];
for qi=1:length(Qs)
  Q = Qs(qi);
  for d=ds
    cost(d,qi) = 0;
    for c=1:maxC
      if length(Qsize(c,d) > 0) % this clique contains Q nodes
	cost(d,qi) = cost(d,qi) + Q^Qsize(c,d)*2^Fsize(c,d);
      end
    end
    %cost_bound(d,qi) = effnumclqs(d) * 8 * Q^(max(Qsize(:,d)));
    cost_bound(d,qi) = (effnumclqs(d)*8) + Q^(max(Qsize(:,d)));
  end
end

qi=2; plot(ds, cost(:,qi), 'o-',  ds, cost_bound(:,qi), '*--');
end


if 0
% convert numbers in cliques into names
for d=1:D
  Fdecode(Fnodes(d)) = d;
end
for c=8:15
  clqs = J.cliques{c};
  fprintf('clique %d: ', c);
  for k=clqs
    if myismember(k, Qnodes)
      fprintf('Q%d ', k)
    elseif myismember(k, Fnodes)
      fprintf('F%d ', Fdecode(k))
    elseif isequal(k, Onode)
      fprintf('O ')
    elseif myismember(k, Qnodes2)
      fprintf('Q%d* ', k-ss)
    else
      error(['unrecognized node ' k])
    end
  end
  fprintf('\n');
end
end
