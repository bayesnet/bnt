function p = adjustable_CPD(CPD)
% ADJUSTABLE_CPD Does this CPD have any adjustable params? (gaussian)
% p = adjustable_CPD(CPD)

p = ~CPD.clamped_mean | ~CPD.clamped_cov | ~CPD.clamped_weights;
