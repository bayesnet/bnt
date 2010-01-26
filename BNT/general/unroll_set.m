function U = unroll_set(S, ss, T)
% UNROLL_SET Make T shifted copies of the set of nodes S in a slice of size ss.
% U = unroll_set(S, ss, T)

offset = repmat(0:ss:(T-1)*ss, [length(S) 1]);
U = repmat(S(:), [1 T]) + offset;
    
