% Example of explaining away from
% http://www.ai.mit.edu/~murphyk/Bayes/bnintro.html#explainaway
%
% Suppose you have to be brainy or smart to get into college.
% B S P(C=1) P(C=2)  1=false 2=true
% 1 1 1.0    0.0 
% 2 1 0.0    1.0
% 1 2 0.0    1.0
% 2 2 0.0    1.0
%
%
% If we observe that you are in college, you must be either brainy or sporty or both.
% If we observre you are in college and sporty, it is less likely you are brainy, 
% since brainy-ness and sporty-ness compete as causal explanations of the effect.

% B  S
%  \/
%   C

B = 1; S = 2; C = 3;
dag = zeros(3,3);
dag([B S], C)=1;
ns = 2*ones(1,3);
bnet = mk_bnet(dag, ns);
bnet.CPD{B} = tabular_CPD(bnet, B, 'CPT', [0.5 0.5]');
bnet.CPD{S} = tabular_CPD(bnet, S, 'CPT', [0.5 0.5]');
CPT = zeros(2,2,2);
CPT(1,1,:) = [1 0];
CPT(2,1,:) = [0 1];
CPT(1,2,:) = [0 1];
CPT(2,2,:) = [0 1];
bnet.CPD{C} = tabular_CPD(bnet, C, 'CPT', CPT);

engine = jtree_inf_engine(bnet);
ev = cell(1,3);
ev{C} = 2;
engine = enter_evidence(engine, ev);
m = marginal_nodes(engine, B);
fprintf('P(B=true|C=true) = %5.3f\n', m.T(2)) % 0.67

ev{S} = 2;
engine = enter_evidence(engine, ev);
m = marginal_nodes(engine, B);
fprintf('P(B=true|C=true,S=true) = %5.3f\n', m.T(2)) % 0.5 = unconditional baseline P(B=true)
