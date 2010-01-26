function pot = enter_cts_evidence_pot(pot, Y, y)
% function pot = enter_cts_evidence_pot(pot, Y, y) cgpot


pot = cg_mom_to_can(pot);
for i=1:pot.dsize
  pot.can{i} = enter_cts_evidence_pot(pot.can{i}, Y, y);
end
