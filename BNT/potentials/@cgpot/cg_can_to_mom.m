function pot = cg_can_to_mom(pot)
% CG_CAN_TO_MOM Convert a CG potential from canonical to moment form, if necessary.
% pot = cg_can_to_mom(pot)

if pot.subtype ~= 'm'
  for i=1:pot.dsize
    pot.mom{i} = cpot_to_mpot(pot.can{i});
  end
  pot.subtype = 'm';
end   
