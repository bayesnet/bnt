function [mpe, ll] = calc_mpe_global(bnet, evidence)
% CALC_MPE_GLOBAL Compute the most probable explanation(s) from the global joint
% [mpe, ll] = calc_mpe_global(bnet, evidence)
%
% mpe(k,i) is the most probable value of node i in the k'th global mode 
% ll is the log likelihood
%
% We assume all nodes are discrete

engine = global_joint_inf_engine(bnet);
engine = enter_evidence(engine, evidence);
S1 = struct(engine); % violate object privacy
S2 = struct(S1.jpot); % joint potential
prob = max(S2.T(:));
modes = find(S2.T(:) == prob);

ens = bnet.node_sizes;
onodes = find(~isemptycell(evidence));
ens(onodes) = 1;
mpe = ind2subv(ens, modes);
for k=1:length(modes)
  for i=onodes(:)'
    mpe(k,i) = evidence{i};
  end
end
ll = log(prob);

mpe = num2cell(mpe);
