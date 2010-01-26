function [engine, clqtoroot] = push(engine, clq, pushdom)
%PUSH_POT push the variables in putshdom which is subset of clq to the clique toword the root and get new engine
%pushdom is pushed variables set
%clq is the index of the clique that pushdom belongs to

clqdom = engine.cliques{clq};
assert( mysubset(pushdom, clqdom));
clqtoroot = parents(engine.jtree, clq);
%sepdom = engine.separator{clq, clqtoroot};
sepdom = engine.separator{clqtoroot, clq};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the strong marginal of the union of pushdom and and the separatordomain and  %
% the corresponding complement                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[margpot, comppot] = complement_pot(engine.clpot{clq}, pushdom);
newsepdom = myunion(pushdom,sepdom);
[margpot,comppot] = complement_pot(engine.clpot{clq}, newsepdom);
engine.clpot{clqtoroot} = direct_combine_pots(engine.clpot{clqtoroot}, margpot);
engine.clpot{clq} = comppot;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of the new separator and separatorpotential of the junction tree %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
engine.seppot{clqtoroot, clq} = direct_combine_pots(engine.seppot{clqtoroot, clq}, margpot);
engine.separator{clqtoroot, clq} = myunion(engine.separator{clqtoroot, clq}, pushdom);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add pushdomain to the clique towards the root %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
engine.cliques{clqtoroot} = myunion(engine.cliques{clqtoroot}, pushdom);

num_cliques = length(engine.cliques);
B = sparse(num_cliques, 1);
for i=1:num_cliques
  B(i, engine.cliques{i}) = 1;
end
engine.cliques_bitv = B;
