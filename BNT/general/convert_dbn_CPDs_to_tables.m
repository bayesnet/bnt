function CPDpot = convert_dbn_CPDs_to_tables(bnet, evidence)
% CONVERT_DBN_CPDS_TO_TABLES Convert CPDs of (possibly instantiated) DBN nodes to tables
% CPDpot = convert_dbn_CPDs_to_tables(bnet, evidence)
%
% CPDpot{n,t} is a table containing P(n,t|pa(n,t), ev)
% All hidden nodes are assumed to be discrete.
% We assume the observed nodes are the same in every slice.
%
% Evaluating the conditional likelihood of long evidence sequences can be very slow,
% so we take pains to vectorize where possible.

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

% special cases: c=child, p=parents, d=discrete, h=hidden, 1sl=1slice
% if c=h=1 then c=d=1, since hidden nodes must be discrete
% c=h c=d p=h p=d 1sl method
% ---------------------------
% 1   1   1   1   -   replicate CPT
% -   1   -   1   -   evaluate CPT on evidence *
% 0   1   1   1   1   dhmm
% 0   0   1   1   1   ghmm
% other               loop
%
% * = any subset of the domain may be observed

% Example where all of the special cases occur - a hierarchical HMM
% where the top layer (G) and leaves (Y) are observed and
% all nodes are discrete except Y.
% (O turns on if Y is an outlier)

% G ---------> G 
% |            |
% v            v
% S  --------> S
% |            |
% v            v
% Y            Y
% ^            ^
% |            |
% O            O

% Evaluating P(yt|St,Ot) is the ghmm case
% Evaluating P(St|S(t-1),gt) is the eval CPT case
% Evaluating P(gt|g(t-1) is the eval CPT case (hdom = [])
% Evaluating P(Ot) is the replicated CPT case

% Cts parents (e.g., inputs) would require an additional special case for speed


  % slices 2..T
  [ss T] = size(evidence);
  self = n+ss;
  ps = parents(bnet.dag, self);
  e = bnet.equiv_class(n, 2);

  if 1
  debug = 0;
  hidden_child = ~obs_bitv(n);
  discrete_child = myismember(n, bnet.dnodes);
  hidden_ps = all(~obs_bitv(ps));
  discrete_ps = mysubset(ps, bnet.dnodes);
  parents_in_same_slice = all(ps > ss);
  
  if hidden_child & discrete_child & hidden_ps & discrete_ps
    CPDpot = helper_repl(bnet, evidence, n, CPDpot, obs_bitv, debug);
  elseif discrete_child & discrete_ps
    CPDpot = helper_eval(bnet, evidence, n, CPDpot, obs_bitv, debug);
  elseif discrete_child & hidden_ps & discrete_ps & parents_in_same_slice
    CPDpot = helper_dhmm(bnet, evidence, n, CPDpot, obs_bitv, debug);
  elseif ~discrete_child & hidden_ps & discrete_ps & parents_in_same_slice
    CPDpot = helper_ghmm(bnet, evidence, n, CPDpot, obs_bitv, debug);
  else
    if debug, fprintf('node %d, slow\n', n); end
    for t=2:T
      CPDpot{n,t} = convert_to_table(bnet.CPD{e}, [ps self], evidence(:,t-1:t));
    end
  end
  end
  
  if 0
  for t=2:T
    CPDpot2{n,t} = convert_to_table(bnet.CPD{e}, [ps self], evidence(:,t-1:t));
    if ~approxeq(CPDpot{n,t}, CPDpot2{n,t})
      fprintf('CPDpot n=%d, t=%d\n',n,t);
      keyboard
    end
  end
  end

  
end




%%%%%%%
function CPDpot = helper_repl(bnet, evidence, n, CPDpot, obs_bitv, debug)

[ss T] = size(evidence);
if debug, fprintf('node %d, repl\n', n); end
e = bnet.equiv_class(n, 2);
CPT = convert_CPD_to_table_hidden_ps(bnet.CPD{e}, []);
CPDpot(n,2:T) = num2cell(repmat(CPT, [1 1 T-1]), [1 2]);



%%%%%%%
function CPDpot = helper_eval(bnet, evidence, n, CPDpot, obs_bitv, debug)

[ss T] = size(evidence);
self = n+ss;
ps = parents(bnet.dag, self);
e = bnet.equiv_class(n, 2);
ns = bnet.node_sizes(:);
% Example: given CPT(p1, p2, p3, p4, c), where p1,p3 are observed
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
if isempty(hdom)
  CPT = CPT(:);
else
  CPT = permute(CPT, [hdom_rel odom_rel]);
  CPT = reshape(CPT, prod(ns(hdom)), prod(ns(odom)));
end
parents_in_same_slice = all(ps > ss);
if parents_in_same_slice
  if debug, fprintf('node %d eval 1 slice\n', n); end
  data = cell2num(evidence(odom-ss,2:T)); %data(i,t) = val of i'th obs parent at t+1
else
  if debug, fprintf('node %d eval 2 slice\n', n); end
  % there's probably a way of vectorizing this...
  data = zeros(length(odom), T-1);
  for t=2:T
    ev = evidence(:,t-1:t);
    ev = ev(:);
    ev2 = ev(odom);
    data(:,t-1) = cat(1, ev2{:});
    %data(:,t-1) = cell2num(ev2);
  end
end
ndx = subv2ind(ns(odom), data'); % ndx(t) encodes data(:,t)
if isempty(hdom)
  CPDpot(n,2:T) = num2cell(CPT(ndx)); % a cell array of floats
else
  CPDpot(n,2:T) = num2cell(CPT(:, ndx), 1); % a cell array of column vectors
end

%%%%%%%
function CPDpot = helper_dhmm(bnet, evidence, n, CPDpot, obs_bitv, debug)

if debug, fprintf('node %d, dhmm\n', n); end
[ss T] = size(evidence);
self = n+ss;
ps = parents(bnet.dag, self);
e = bnet.equiv_class(n, 2);
ns = bnet.node_sizes(:);
CPT = CPD_to_CPT(bnet.CPD{e});
CPT = reshape(CPT, [prod(ns(ps)) ns(self)]); % what if no parents?
%obslik = mk_dhmm_obs_lik(cell2num(evidence(n,2:T)), CPT);
obslik = eval_pdf_cond_multinomial(cell2num(evidence(n,2:T)), CPT);
CPDpot(n,2:T) = num2cell(obslik, 1);


%%%%%%%
function CPDpot = helper_ghmm(bnet, evidence, n, CPDpot, obs_bitv, debug)

if debug, fprintf('node %d, ghmm\n', n); end
[ss T] = size(evidence);
e = bnet.equiv_class(n, 2);
S = struct(bnet.CPD{e}); 
ev2 = cell2num(evidence(n,2:T));
%obslik = mk_ghmm_obs_lik(ev2, S.mean, S.cov);
obslik = eval_pdf_cond_gauss(ev2, S.mean, S.cov);
CPDpot(n,2:T) = num2cell(obslik, 1);

