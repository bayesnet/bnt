function CPD = generic_CPD(clamped)
% GENERIC_CPD Virtual constructor for generic CPD
% CPD = discrete_CPD(clamped)

if nargin < 1, clamped = 0; end

CPD.clamped = clamped;
CPD = class(CPD, 'generic_CPD');
