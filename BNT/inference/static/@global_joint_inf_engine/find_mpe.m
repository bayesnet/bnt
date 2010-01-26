function [mpe, ll] = find_mpe(engine, evidence)
% FIND_MPE_GLOBAL Compute the most probable explanation(s) from the global joint
% [mpe, ll] = find_mpe(engine, evidence)
%
% mpe(k,i) is the most probable value of node i in the k'th global mode  (cell array)
%
% We assume all nodes are discrete

%engine = global_joint_inf_engine(bnet);
bnet = bnet_from_engine(engine);
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
