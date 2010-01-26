function CPD = discrete_CPD(clamped, dom_sizes)
% DISCRETE_CPD Virtual constructor for generic discrete CPD
% CPD = discrete_CPD(clamped, dom_sizes)

CPD.dom_sizes = dom_sizes;
CPD = class(CPD, 'discrete_CPD', generic_CPD(clamped));
