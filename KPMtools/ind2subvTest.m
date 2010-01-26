function sub = ind2subvTest(siz,index)
%IND2SUBV   Subscript vector from linear index.
% IND2SUBV(SIZ,IND) returns a vector of the equivalent subscript values 
% corresponding to a single index into an array of size SIZ.
% If IND is a vector, then the result is a matrix, with subscript vectors
% as rows.

sub = ind2subvMinka(siz, index);
subKPM = ind2subvKPM(siz, index);
assert(isequal(sub, subKPM))
