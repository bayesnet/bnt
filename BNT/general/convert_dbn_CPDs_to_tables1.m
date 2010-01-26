function CPDpot = convert_dbn_CPDs_to_tables1(bnet, evidence)
% CONVERT_DBN_CPDS_TO_TABLES Convert CPDs of (possibly instantiated) DBN nodes to tables
% CPDpot = convert_dbn_CPDs_to_tables(bnet, evidence)
%
% CPDpot{n,t} is a table containing P(n,t|pa(n,t), ev)
% All hidden nodes are assumed to be discrete
% We assume the observed nodes are the same in every slice
%
% Evaluating the conditional likelihood of the evidence can be very slow,
% so we take pains to vectorize where possible, i.e., we try to avoid
% calling convert_to_table

[ss T] = size(evidence);
%obs_bitv = ~isemptycell(evidence(:));
obs_bitv = zeros(1, 2*ss);
obs_bitv(bnet.observed) = 1;
obs_bitv(bnet.observed+ss) = 1;

ns = bnet.node_sizes(:);
CPDpot = cell(ss,T); 

for n=1:ss
  % slice 1
  t = 1;
  ps = parents(bnet.dag, n);
  e = bnet.equiv_class(n, 1);
  if ~any(obs_bitv(ps))
    CPDpot{n,t} = convert_CPD_to_table_hidden_ps(bnet.CPD{e}, evidence{n,t});
  else
    CPDpot{n,t} = convert_to_table(bnet.CPD{e}, [ps n], evidence(:,1));
  end
  
  % slices 2..T
  debug = 1;
  if ~obs_bitv(n)
    CPDpot = helper_hidden_child(bnet, evidence, n, CPDpot, obs_bitv, debug);
  else
    CPDpot = helper_obs_child(bnet, evidence, n, CPDpot, obs_bitv, debug);
  end
end

if 0
CPDpot2 = convert_dbn_CPDs_to_tables_slow(bnet, evidence);
for t=1:T
  for n=1:ss
    if ~approxeq(CPDpot{n,t}, CPDpot2{n,t})
      fprintf('CPDpot n=%d, t=%d\n',n,t);
      keyboard
    end
  end
end
end


% special cases: c=child, p=parents, d=discrete, h=hidden, 1=1slice
% if c=h=1 then c=d=1, since hidden nodes must be discrete
% c=h c=d p=h p=d p=1 method
% ---------------------------
% 1   1   1   1   -   replicate CPT
% 0   1   1   1   1   dhmm
% 0   0   1   1   1   ghmm
% -   1   -   1   -   evaluate CPT on evidence
% other               loop

%%%%%%%
function CPDpot = helper_hidden_child(bnet, evidence, n, CPDpot, obs_bitv, debug)

[ss T] = size(evidence);
self = n+ss;
ps = parents(bnet.dag, self);
e = bnet.equiv_class(n, 2);
ns = bnet.node_sizes(:);
if ~any(obs_bitv(ps)) % all parents are hidden (hence discrete)
  if debug, fprintf('node %d is hidden, all ps are hidden\n', n); end
  if myismember(n, bnet.dnodes) 
    %CPT = CPD_to_CPT(bnet.CPD{e});
    %CPT = reshape(CPT, [prod(ns(ps)) ns(self)]);
    CPT = convert_CPD_to_table_hidden_ps(bnet.CPD{e}, []);
    CPDpot(n,2:T) = num2cell(repmat(CPT, [1 1 T-1]), [1 2]);
  else
    error(['hidden cts node disallowed'])
  end
else % some parents are observed - slow
  if mysubset(ps, bnet.dnodes) % all parents are discrete
    % given CPT(p1, p2, p3, p4, c), where p1,p3 are observed
    % we create CPT([p2 p4 c], [p1 p3]).
    % We then convert all observed p1,p3 into indices ndx
    % and return CPT(:, ndx)
    CPT = CPD_to_CPT(bnet.CPD{e});
    domain = [ps self];
    % if dom is [3 7 8] and 3,8 are observed, odom_rel = [1 3], hdom_rel = 2,
    % odom = [3 8], hdom = 7
    odom_rel = find(obs_bitv(domain));
    hdom_rel = find(~obs_bitv(domain));
    odom = domain(odom_rel);
    hdom = domain(hdom_rel);
    CPT = permute(CPT, [hdom_rel odom_rel]);
    CPT = reshape(CPT, prod(ns(hdom)), prod(ns(odom)));
    parents_in_same_slice = all(ps > ss);
    if parents_in_same_slice
      if debug, fprintf('node %d is hidden, some ps are obs, all ps discrete, 1 slice\n', n); end
      data = cell2num(evidence(odom-ss,2:T)); %data(i,t) = val of i'th obs parent at t+1
    else
      if debug, fprintf('node %d is hidden, some ps are obs, all ps discrete, 2 slice\n', n); end
      data = zeros(length(odom), T-1);
      for t=2:T
	ev = evidence(:,t-1:t);
	data(:,t-1) = cell2num(ev(odom));
      end
    end
    ndx = subv2ind(ns(odom), data'); % ndx(t) encodes data(:,t)
    CPDpot(n,2:T) = num2cell(CPT(:, ndx), [1 2]);
  else % some parents are cts - v slow
    if debug, fprintf('node %d is hidden, some ps are obs, some ps cts\n', n); end
    for t=2:T
      CPDpot{n,t} = convert_to_table(bnet.CPD{e}, [ps self], evidence(:,t-1:t));
    end
  end
end
  
%%%%%%%
function CPDpot = helper_obs_child(bnet, evidence, n, CPDpot, obs_bitv, debug)

[ss T] = size(evidence);
self = n+ss;
ps = parents(bnet.dag, self);
e = bnet.equiv_class(n, 2);
ns = bnet.node_sizes(:);
if ~any(obs_bitv(ps)) % all parents are hidden
  parents_in_same_slice = all(ps > ss);
  if parents_in_same_slice
    if debug, fprintf('node %d is obs, all ps are hidden, 1 slice\n', n); end
    ps1 = ps - ss;
    if myismember(n, bnet.dnodes) 
      CPT = CPD_to_CPT(bnet.CPD{e});
      CPT = reshape(CPT, [prod(ns(ps)) ns(self)]); % what if no parents?
      obslik = eval_pdf_cond_multinomial(cell2num(evidence(n,2:T)), CPT);
      CPDpot(n,2:T) = num2cell(obslik, 1);
    else
      S = struct(bnet.CPD{e}); 
      obslik = eval_pdf_cond_gauss(cell2num(evidence(n,2:T)), S.mean, S.cov);
      CPDpot(n,2:T) = num2cell(obslik, 1);
    end
  else % parents span 2 slices - slow
    if debug, fprintf('node %d is obs, all ps are hidden , 2 slice\n', n); end
    for t=2:T
      CPDpot{n,t} = convert_to_table(bnet.CPD{e}, [ps self], evidence(:,t-1:t));
    end
  end
else 
  if isempty(ps) % observed root
    if debug, fprintf('node %d is obs, no ps\n', n); end
    CPT = CPD_to_CPT(bnet.CPD{e});
    data = cell2num(evidence(n,2:T));
    CPDpot(n,2:T) = CPT(data);
  else  % some parents are observed  - slow
    if debug, fprintf('node %d is obs, some ps are obs\n', n); end
    for t=2:T
      CPDpot{n,t} = convert_to_table(bnet.CPD{e}, [ps self], evidence(:,t-1:t));
    end
  end
end
