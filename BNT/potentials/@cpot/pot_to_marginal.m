function m = pot_to_marginal(pot)
% POT_TO_MARGINAL Convert a cpot to a marginal structure.
% m = pot_to_marginal(pot)

mom = cpot_to_mpot(pot);
m = pot_to_marginal(mom);
