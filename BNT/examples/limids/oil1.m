% oil wildcatter influence diagram in Cowell et al p172

% T = test for oil?
% UT = utility (negative cost) of testing
% O = amount of oil = Dry, Wet or Soaking
% R = results of test = NoStrucure, OpenStructure, ClosedStructure or NoResult
% D = drill?
% UD = utility of drilling

% Decision sequence = T R D O

T = 1; UT = 2; O = 3; R = 4; D = 5; UD = 6;
N = 6;
dag = zeros(N);
dag(T, [UT R D]) = 1;
dag(O, [R UD]) = 1;
dag(R, D) = 1;
dag(D, UD) = 1;

ns = zeros(1,N);
ns(O) = 3; ns(R) = 4; ns(T) = 2; ns(D) = 2; ns(UT) = 1; ns(UD) = 1;

limid = mk_limid(dag, ns, 'chance', [O R], 'decision', [T D], 'utility', [UT UD]);

limid.CPD{O} = tabular_CPD(limid, O, [0.5 0.3 0.2]);
tbl = [0.6 0 0.3 0 0.1 0  0.3 0 0.4 0 0.4 0  0.1 0 0.3 0 0.5 0  0 1 0 1 0 1];
limid.CPD{R} = tabular_CPD(limid, R, tbl);

limid.CPD{UT} = tabular_utility_node(limid, UT,  [-10 0]);
limid.CPD{UD} = tabular_utility_node(limid, UD, [-70 50 200  0 0 0]);

if 1
  % start with uniform policies
  limid.CPD{T} = tabular_decision_node(limid, T);
  limid.CPD{D} = tabular_decision_node(limid, D);
else
  % hard code optimal policies
  limid.CPD{T} = tabular_decision_node(limid, T, [1.0 0.0]);        
  a = 0.5; b = 1-a; % arbitrary value
  tbl = myreshape([0 a 1 a 1 a a a  1 b 0 b 0 b b b], ns([T R D]));
  limid.CPD{D} = tabular_decision_node(limid, D,  tbl);
end

%fname = '/home/cs/murphyk/matlab/Misc/loopybel.txt';

engines = {};
engines{end+1} = global_joint_inf_engine(limid);
engines{end+1} = jtree_limid_inf_engine(limid);
%engines{end+1} = belprop_inf_engine(limid, 'max_iter', 3*N, 'filename', fname);

exact = [1 2];
%approx = 3;
approx = [];

E = length(engines);
strategy = cell(1, E);
MEU = zeros(1, E);
for e=1:E
  [strategy{e}, MEU(e)] = solve_limid(engines{e});
  MEU
end
MEU

for e=exact(:)'
  assert(approxeq(MEU(e), 22.5))
  % U(T=yes)  U(T=no)
  % 1         0
  assert(argmax(strategy{e}{T}) == 1); % test = yes
  t = 1; % test = yes
  % strategy{D} T       R      U(D=yes=1)  U(D=no=2)
  %             1=yes   1=noS  0           1         Don't drill
  %             2=no    1=noS  1           0
  %             1=yes   2=opS  1           0
  %             2=no    2=opS  1           0
  %             1=yes   3=clS  1           0
  %             2=no    3=clS  1           0
  %             1=yes   4=unk  1           0
  %             2=no    4=unk  1           0
  
  for r=[2 3] % OpS, ClS
    assert(argmax(squeeze(strategy{e}{D}(t,r,:))) == 1); % drill = yes
  end
  r = 1; % noS
  assert(argmax(squeeze(strategy{e}{D}(t,r,:))) == 2); % drill = no
end


for e=approx(:)'
  approxeq(strategy{exact(1)}{T}, strategy{e}{T})
  approxeq(strategy{exact(1)}{D}, strategy{e}{D})
end
