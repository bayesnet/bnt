function pot = rescale_pot(pot, s)
% RESCALE_POT Add a constant to the cpot scale factor.
% pot = rescale_pot(pot, s)

pot.g = pot.g + s;
