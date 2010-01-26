function pot = rescale_pot(pot, s)
% RESCALE_POT Add a constant to the mpot scale factor.
% pot = rescale_pot(pot, s)

pot.p = pot.p*s;
