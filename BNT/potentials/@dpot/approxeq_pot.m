function p = approxeq_pot(A, B, tol)

if nargin < 3, tol = 1e-3; end

p = approxeq(A.T, B.T, tol);
