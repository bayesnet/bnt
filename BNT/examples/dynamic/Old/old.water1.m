% Compare the speeds of various inference engines on the water DBN

[bnet, onodes] = mk_water_dbn;

T = 3;

engine = {};
engine{end+1} = jtree_unrolled_dbn_inf_engine(bnet, T, onodes);
engine{end+1} = hmm_inf_engine(bnet, onodes);
engine{end+1} = frontier_inf_engine(bnet, onodes);
engine{end+1} = jtree_dbn_inf_engine(bnet, onodes);
engine{end+1} = bk_inf_engine(bnet, 'exact', onodes);

engine{end+1} = bk_inf_engine(bnet, 'ff', onodes);
engine{end+1} = bk_inf_engine(bnet, { [1 2], [3 4 5 6], [7 8] }, onodes);

N = length(engine);
exact = 1:5;


filter = 0;
err = cmp_inference(bnet, onodes, engine, exact, T, filter);

% elapsed times for enter_evidence  (matlab 5.3 on PIII with 256MB running Redhat linux)  

% T = 5, 4/20/00
%     0.6266 unrolled *
%     0.3490 hmm *
%     1.1743 frontier
%     1.4621 old frontier
%     0.3270 fast frontier *
%     1.3926 jtree 
%     1.3790 bk 
%     0.4916 fast bk
%     0.4190 fast bk compiled
%     0.3574 fast jtree *


err = cmp_learning(bnet, onodes, engine, exact, T);

% elapsed times for learn_params_dbn_em (matlab 5.3 on PIII with 256MB running Redhat linux)  

% T = 5, 2cases, 2 iter, 4/20/00
% 3.5750 unrolled
% 3.7475 hmm
% 2.1452 fast frontier
% 2.5724 fast bk compiled
% 2.3387 fast jtree
