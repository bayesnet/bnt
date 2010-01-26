seed = 0;
rand('state', seed);
randn('state', seed);

nrows = 3;
ncols = 3;
npixels = nrows*ncols;

% we number pixels in transposed raster scan order (top to bottom, left to right)

% hidden var
HV = reshape(1:npixels, nrows, ncols);
% observed var
OV = reshape(1:npixels, nrows, ncols) + length(HV(:));

% observed factor
OF = reshape(1:npixels, nrows, ncols);
% vertical edge factor VEF(i,j) is the factor for edge HV(i,j) - HV(i+1,j)
VEF = reshape((1:(nrows-1)*ncols), nrows-1, ncols) + length(OF(:));
% horizontal edge factor HEF(i,j) is the factor for edge HV(i,j) - HV(i,j+1)
HEF = reshape((1:nrows*(ncols-1)), nrows, ncols-1) + length(OF(:)) + length(VEF(:));

nvars = length(HV(:))+length(OV(:));
assert(nvars == 2*npixels);
nfac = length(OF(:)) + length(VEF(:)) + length(HEF(:));

K = 2; % number of discrete values for the hidden vars
%O = 1; % each observed pixel is a scalar
O = 2; % each observed pixel is binary

factors = cell(1,3);

% hidden states generate observed 0 or 1 plus noise
%factors{2} = cond_gauss1_kernel(K, O, 'mean', [0 1], 'cov', [0.1 0.1]);
pnoise = 0.2;
factors{1} = tabular_kernel([K O], [1-pnoise pnoise; pnoise 1-pnoise]);
ofactor = 1;

% encourage compatibility between neighboring vertical pixels
factors{2} = tabular_kernel([K K], [0.8 0.2; 0.2 0.8]);
vedge_factor = 2;

%% no constraint between neighboring horizontal pixels
%factors{3} = tabular_kernel([K K], [0.5 0.5; 0.5 0.5]);

factors{3} = tabular_kernel([K K], [0.8 0.2; 0.2 0.8]);
hedge_factor = 3;



factor_ndx = zeros(1, 3);
G = zeros(nvars, nfac);
ns = [K*ones(1,length(HV(:))) O*ones(1,length(OV(:)))];

N = length(ns);
%cnodes = OV(:);
cnodes = [];
dnodes = 1:N;

for i=1:nrows
  for j=1:ncols
    G([HV(i,j), OV(i,j)], OF(i,j)) = 1;
    factor_ndx(OF(i,j)) = ofactor;

    if i < nrows
      G(HV(i:i+1,j), VEF(i,j)) = 1;
      factor_ndx(VEF(i,j)) = vedge_factor;
    end

    if j < ncols
      G(HV(i,j:j+1), HEF(i,j)) = 1;
      factor_ndx(HEF(i,j)) = hedge_factor;
    end

  end
end


fg = mk_fgraph(G, ns, factors, 'discrete', dnodes, 'equiv_class', factor_ndx);

if 1
  % make image with vertical stripes
  I = zeros(nrows, ncols);
  for j=1:2:ncols
    I(:,j) = 1;
  end
else
  % make image with square in middle
  I = zeros(nrows, ncols);
  I(3:6,3:6) = 1;
end

  
% corrupt image
O = mod(I + (rand(nrows,ncols)> (1-pnoise)), 2);

maximize = 1;
engine = belprop_fg_inf_engine(fg, 'maximize', maximize, 'max_iter', npixels*5);

evidence = cell(1, nvars);
onodes = OV(:);
evidence(onodes) = num2cell(O+1); % values must be in range {1,2}

engine = enter_evidence(engine, evidence);

for i=1:nrows
  for j=1:ncols
    m = marginal_nodes(engine, HV(i,j));
    Ihat(i,j) = argmax(m.T)-1;
  end
end

Ihat
