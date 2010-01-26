function s = checkpsd(s)

if (any(isnan(s) | isinf(s) | ~isreal(s)))
  warning('S contains complex numbers, Inf, or NaN'); 
end
% Drop any negative eigenvalues.
[V, D] = eig(full(s));
d = real(diag(D));
if (any(d < 0))
  warning(sprintf(['S is not positive semidefinite (min. eig. =' ...
		   ' %0.5g); projecting.'], min(d)));
  d(find(d < 0)) = 0;
  D = diag(d);
  s = V * D * V';
end
