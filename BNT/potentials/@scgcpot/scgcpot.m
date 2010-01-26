function pot = scgcpot(cheadsize, ctailsize, p, A, B, C)
% SCGCPOT Make a base object of stable conditional gaussian potential.
% pot = scgcpot(cheadsize, ctailsize, p, A, B, C)
%
% cheadsize is the demension of head nodes.
% ctailsize is the demension of tail nodes.
% r = cheadsize, s = ctailsize
% p is discrete probability.
% A is table of r*1 vectors;
% B is r*s matrices
% C is r*r positive semidefinite symmetric matrices

if nargin < 3
    p = 1; 
end
if nargin < 4
    A = zeros(cheadsize,1); 
end
if nargin < 5
    B = zeros(cheadsize,ctailsize); 
end
if nargin < 6
    C = zeros(cheadsize,cheadsize); 
end

if isempty(A)
    A = zeros(cheadsize,1); 
end
if isempty(B)
    B = zeros(cheadsize,ctailsize); 
end
if isempty(C)
    C = zeros(cheadsize,cheadsize); 
end
  
pot.cheadsize = cheadsize;
pot.ctailsize = ctailsize;

pot.p = p;
pot.A = A;
pot.B = B;
pot.C = C;
%if cheadsize == 0
%   pot.A = [];
%end
%if ctailsize == 0
%    pot.B = [];
%end
pot = class(pot, 'scgcpot');
