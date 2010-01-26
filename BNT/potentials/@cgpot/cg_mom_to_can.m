function pot = cg_mom_to_can(pot)
% CG_MOM_TO_CAN Convert a CG potential from moment to canonical form, if necessary.
% pot = cg_mom_to_can(pot)

if pot.subtype ~= 'c'
  for i=1:pot.dsize
    pot.can{i} = mpot_to_cpot(pot.mom{i});
  end
  pot.subtype = 'c';
end   
