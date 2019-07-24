function ndx = subv2indTest(siz,sub)
%SUBV2IND   Linear index from subscript vector.
% SUBV2IND(SIZ,SUB) returns an equivalent single index corresponding to a
% subscript vector for an array of size SIZ.
% If SUB is a matrix, with subscript vectors as rows, then the result is a 
% column vector.
%
% This is the opposite of IND2SUBV, so that
%   SUBV2IND(SIZ,IND2SUBV(SIZ,IND)) == IND.
%
% See also IND2SUBV, SUB2IND.

ndx = subv2indMinka(siz,sub);
ndx2 = subv2indKPM(siz,sub);
assert(isequalKPM(ndx,ndx2))
