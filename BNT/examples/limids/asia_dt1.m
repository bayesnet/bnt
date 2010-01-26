% decision theoretic version of asia network
% Cowell et al, p177
% We explicitely add the no-forgetting arcs.

Smoking = 1;
VisitToAsia = 2;
Bronchitis = 3;
LungCancer = 4;
TB = 5;
Do_Xray = 6;
TBorCancer = 7;
Util_Xray = 8;
Dys = 9;
posXray = 10;
Do_Hosp = 11;
Util_Hosp = 12;

n = 12;
dag = zeros(n);
dag(Smoking, [Bronchitis LungCancer]) = 1;
dag(VisitToAsia, [TB Do_Xray Do_Hosp]) = 1;
dag(Bronchitis, Dys) = 1;
dag(LungCancer, [Util_Hosp TBorCancer]) = 1;
dag(TB, [Util_Hosp TBorCancer Util_Xray]) = 1;
dag(Do_Xray, [posXray Util_Xray Do_Hosp]) = 1;
dag(TBorCancer, [Dys posXray]) = 1;
dag(Dys, Do_Hosp) = 1;
dag(posXray, Do_Hosp) = 1;
dag(Do_Hosp, Util_Hosp) = 1;

dnodes = [Do_Xray Do_Hosp];
unodes = [Util_Xray Util_Hosp];
cnodes = mysetdiff(1:n, [dnodes unodes]); % chance nodes
ns = 2*ones(1,n);
ns(unodes) = 1;
limid = mk_limid(dag, ns, 'chance', cnodes, 'decision', dnodes, 'utility', unodes);

% 1 = yes, 2 = no
limid.CPD{VisitToAsia} = tabular_CPD(limid, VisitToAsia, [0.01 0.99]);
limid.CPD{Bronchitis} = tabular_CPD(limid, Bronchitis, [0.6 0.3  0.4 0.7]);
limid.CPD{Dys} = tabular_CPD(limid, Dys, [0.9 0.7 0.8 0.1  0.1 0.3 0.2 0.9]);
limid.CPD{TBorCancer} = tabular_CPD(limid, TBorCancer, [1 1 1 0  0 0 0 1]);

limid.CPD{LungCancer} = tabular_CPD(limid, LungCancer, [0.1 0.01  0.9 0.99]);
limid.CPD{Smoking} = tabular_CPD(limid, Smoking, [0.5 0.5]);
limid.CPD{TB} = tabular_CPD(limid, TB, [0.05 0.01  0.95 0.99]);
limid.CPD{posXray} = tabular_CPD(limid, posXray, [0.98 0.5 0.05 0.5  0.02 0.5 0.95 0.5]);

limid.CPD{Util_Hosp} = tabular_utility_node(limid, Util_Hosp, [180 120 160 15  2 4 0 40]);
limid.CPD{Util_Xray} = tabular_utility_node(limid, Util_Xray, [0 1 10 10]);

for i=dnodes(:)'
  limid.CPD{i} = tabular_decision_node(limid, i);
end

engines = {};
engines{end+1} = global_joint_inf_engine(limid);
engines{end+1} = jtree_limid_inf_engine(limid);
%engines{end+1} = belprop_inf_engine(limid);

exact = [1 2];
%approx = 3;
approx = [];


NE = length(engines);
MEU = zeros(1, NE);
niter = zeros(1, NE);
strategy = cell(1, NE);

tol = 1e-2;
for e=1:length(engines)
  [strategy{e}, MEU(e), niter(e)] = solve_limid(engines{e});
end

for e=exact(:)'
  assert(approxeq(MEU(e), 47.49, tol))
  assert(isequal(strategy{e}{Do_Xray}(:)', [1 0 0 1]))
  
  % Check the hosptialize strategy is correct (p180)
  % We assume the patient has not been to Asia and therefore did not have an Xray.
  % In this case it is optimal not to hospitalize regardless of whether the patient has
  % dyspnoea or not (and of course regardless of the value of pos_xray).
  asia = 2;
  do_xray = 2;
  for dys = 1:2
    for pos_xray = 1:2
      assert(argmax(squeeze(strategy{e}{Do_Hosp}(asia, do_xray, dys, pos_xray, :))) == 2)
    end
  end
end


for e=approx(:)'
  approxeq(strategy{exact(1)}{Do_Xray}, strategy{e}{Do_Xray})
  approxeq(strategy{exact(1)}{Do_Hosp}, strategy{e}{Do_Hosp})
end

 

