function p = approxeq_pot(A, B, tol)

if nargin < 3, tol = 1e-3; end

p = approxeq(A.p, B.p, tol) & approxeq(A.u, B.u, tol);
