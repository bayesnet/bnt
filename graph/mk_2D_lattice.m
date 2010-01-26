function G = mk_2D_lattice(nrows, ncols, con)
% MK_2D_LATTICE Return adjacency matrix for nearest neighbor connected 2D lattice
% G = mk_2D_lattice(nrows, ncols, con)
% G(k1, k2) = 1 iff k1=(i1,j1) is a neighbor of k2=(i2,j2) 
% (Two pixels are neighbors if their Euclidean distance is less than r.)
% Default connectivity = 4.
%
% WE ASSUME NO WRAP AROUND. 
%
% This is the neighborhood as a function of con:
%
%  con=4,r=1  con=8,r=sqrt(2)   con=12,r=2   con=24,r=sqrt(8)
%  nn            2nd order        4th order
%                                  x         x x x x x
%    x          x x x            x x x       x x x x x
%  x o x        x o x          x x o x x     x x o x x
%    x          x x x            x x x       x x x x x
%                                  x         x x x x x
%
% Examples:
% Consider a 3x4 grid
%  1  4  7  10
%  2  5  8  11
%  3  6  9  12
%  
% 4-connected:
% G=mk_2D_lattice(3,4,4);
% find(G(1,:)) = [2 4]
% find(G(5,:)) = [2 4 6 8]
%
% 8-connected:
% G=mk_2D_lattice(3,4,8);
% find(G(1,:)) = [2 4 5]
% find(G(5,:)) = [1 2 3 4 6 7 8 9]

% meshgrid trick due to Temu Gautama (temu@neuro.kuleuven.ac.be)

if nargin < 3, con = 4; end

switch con,
 case 4, r = 1;
 case 8, r = sqrt(2);
 case 12, r = 2;
 case 24, r = sqrt(8);
 otherwise, error(['unrecognized connectivity ' num2str(con)])
end


npixels = nrows*ncols;

[x y]=meshgrid(1:ncols, 1:nrows);
M = [x(:) y(:)];
M1 = repmat(reshape(M',[1 2 npixels]),[npixels 1 1]);
M2 = repmat(M,[1 1 npixels]);
%D = squeeze(sum(abs(M1-M2),2)); % Manhattan distance
M3 = M1-M2;
D = sqrt(squeeze(M3(:,1,:)) .^2 + squeeze(M3(:,2,:)) .^2); % Euclidean distance
G = reshape(D <= r,npixels,npixels);
G = setdiag(G, 0);
