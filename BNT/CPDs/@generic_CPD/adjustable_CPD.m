function p = adjustable_CPD(CPD)
% ADJUSTABLE_CPD Does this CPD have any adjustable params? (generic)
% p = adjustable_CPD(CPD)
   
p = ~CPD.clamped;
