function G = mk_2D_lattice_slow(nrows, ncols, wrap_around) 
% MK_2D_LATTICE Return adjacency matrix for 4-nearest neighbor connected 2D lattice
% G = mk_2D_lattice(nrows, ncols, wrap_around)
% G(k1, k2) = 1 iff k1=(i1,j1) is connected to k2=(i2,j2)
%
% If wrap_around = 1, we use toroidal boundary conditions (default = 0)
%
% Nodes are assumed numbered as in the following 3x3 lattice
%   1 4 7
%   2 5 8 
%   3 6 9
%
% e.g., G = mk_2D_lattice(3, 3, 0) returns
%   0 1 0 1 0 0 0 0 0 
%   1 0 1 0 1 0 0 0 0 
%   0 1 0 0 0 1 0 0 0 
%   1 0 0 0 1 0 1 0 0 
%   0 1 0 1 0 1 0 1 0 
%   0 0 1 0 1 0 0 0 1 
%   0 0 0 1 0 0 0 1 0 
%   0 0 0 0 1 0 1 0 1 
%   0 0 0 0 0 1 0 1 0 
% so find(G(5,:)) = [2 4 6 8] 
% but find(G(1,:)) = [2 4]
%
% Using wrap around, G = mk_2D_lattice(3, 3, 1), we get
%   0 1 1 1 0 0 1 0 0 
%   1 0 1 0 1 0 0 1 0 
%   1 1 0 0 0 1 0 0 1 
%   1 0 0 0 1 1 1 0 0 
%   0 1 0 1 0 1 0 1 0 
%   0 0 1 1 1 0 0 0 1 
%   1 0 0 1 0 0 0 1 1 
%   0 1 0 0 1 0 1 0 1 
%   0 0 1 0 0 1 1 1 0 
% so find(G(5,:)) = [2 4 6 8] 
% and find(G(1,:)) = [2 3 4 7]

if nargin < 3, wrap_around = 0; end

% M contains the number of each cell e.g.
%   1 4 7
%   2 5 8 
%   3 6 9
% North neighbors (assuming wrap around) are
%   3 6 9
%   1 4 7
%   2 5 8
% Without wrap around, they are
%   1 4 7
%   1 4 7
%   2 5 8
% The first row is arbitrary, since pixels at the top have no north neighbor.

if nrows==1
  G = zeros(1, ncols);
  for i=1:ncols-1
    G(i,i+1) = 1;
    G(i+1,i) = 1;
  end
  if wrap_around
    G(1,ncols) = 1;
    G(ncols,1) = 1;
  end
  return;
end

  
npixels = nrows*ncols;

N = 1; E = 2; S = 3; W = 4;
if wrap_around
  rows{N} = [nrows 1:nrows-1]; cols{N} = 1:ncols;
  rows{E} = 1:nrows; cols{E} = [2:ncols 1];
  rows{S} = [2:nrows 1]; cols{S} = 1:ncols;
  rows{W} = 1:nrows; cols{W} = [ncols 1:ncols-1];
else
  rows{N} = [1 1:nrows-1]; cols{N} = 1:ncols;
  rows{E} = 1:nrows; cols{E} = [1 1:ncols-1];
  rows{S} = [2:nrows nrows]; cols{S} = 1:ncols;
  rows{W} = 1:nrows; cols{W} = [2:ncols ncols];
end

M = reshape(1:npixels, [nrows ncols]);
nbrs = cell(1, 4);
for i=1:4
  nbrs{i} = M(rows{i}, cols{i});
end


G = zeros(npixels, npixels);
if wrap_around
  for i=1:4
    if 0
      % naive
      for p=1:npixels
	G(p, nbrs{i}(p)) = 1;
      end
    else
      % vectorized
      ndx2 = sub2ind([npixels npixels], 1:npixels, nbrs{i}(:)');
      G(ndx2) = 1;
    end
  end
else
  i = N;
  mask = ones(nrows, ncols);
  mask(1,:) = 0; % pixels in row 1 have no nbr to the north
  ndx = find(mask);
  ndx2 = sub2ind([npixels npixels], ndx, nbrs{i}(ndx));
  G(ndx2) = 1;

  i = E;
  mask = ones(nrows, ncols);
  mask(:,ncols) = 0;
  ndx = find(mask);
  ndx2 = sub2ind([npixels npixels], ndx, nbrs{i}(ndx));
  G(ndx2) = 1;

  i = S;
  mask = ones(nrows, ncols);
  mask(nrows,:)=0;
  ndx = find(mask);
  ndx2 = sub2ind([npixels npixels], ndx, nbrs{i}(ndx));
  G(ndx2) = 1;
  
  i = W;
  mask = ones(nrows, ncols);
  mask(:,1)=0;
  ndx = find(mask);
  ndx2 = sub2ind([npixels npixels], ndx, nbrs{i}(ndx));
  G(ndx2) = 1;
end

G = setdiag(G, 0);
