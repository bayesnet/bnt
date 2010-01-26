seed = 0;
rand('state', seed);
randn('state', seed);

nrows = 5;
ncols = 5;
npixels = nrows*ncols;

% we number pixels in transposed raster scan order (top to bottom, left to right)

% H(i,j) is the number of the hidden node at (i,j)
H = reshape(1:npixels, nrows, ncols);
% O(i,j) is the number of the obsevred node at (i,j)
O = reshape(1:npixels, nrows, ncols) + length(H(:));


% Make a Bayes net where each hidden pixel generates an observed pixel
% but there are no connections between the hidden pixels.
% We use this just to generate noisy versions of known images.
N = 2*npixels;
dag = zeros(N);
for i=1:nrows
  for j=1:ncols
    dag(H(i,j), O(i,j)) = 1;
  end
end


K = 2; % number of discrete values for the hidden vars
ns = ones(N,1);
ns(H(:)) = K;
ns(O(:)) = 1;


% make image with vertical stripes
I = zeros(nrows, ncols);
for j=1:2:ncols
  I(:,j) = 1;
end

% each "hidden" node will be instantiated to the pixel in the known image
% each observed node has conditional Gaussian distribution
eclass = ones(1,N);
%eclass(H(:)) = 1;
%eclass(O(:)) = 2;
eclass(H(:)) = 1:npixels;
eclass(O(:)) = npixels+1;
bnet = mk_bnet(dag, ns, 'discrete', H(:), 'equiv_class', eclass);


%bnet.CPD{1} = tabular_CPD(bnet, H(1), 'CPT', normalise(ones(1,K)));
for i=1:nrows
  for j=1:ncols
    bnet.CPD{H(i,j)} = root_CPD(bnet, H(i,j), I(i,j) + 1);
  end
end

% If H(i,j)=1, O(i,j)=+1 plus noise
% If H(i,j)=2, O(i,j)=-1 plus noise
sigma = 0.5;
bnet.CPD{eclass(O(1,1))} = gaussian_CPD(bnet, O(1,1), 'mean', [1 -1], 'cov', reshape(sigma*ones(1,K), [1 1 K]));
ofactor = bnet.CPD{eclass(O(1,1))};
%ofactor = gaussian_CPD('self', 2, 'dps', 1, 'cps', [], 'sz', [K O], 'mean', [1 -1], 'cov', reshape(sigma*ones(1,K), [1 1 K)));


data = sample_bnet(bnet);
img = reshape(data(O(:)), nrows, ncols)




%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Now create MRF represented as a factor graph to try and recover the scene

% VEF(i,j) is the number of the factor for the vertical edge between HV(i,j) - HV(i+1,j)
VEF = reshape((1:(nrows-1)*ncols), nrows-1, ncols);
% HEF(i,j) is the number of the factor for the horizontal edge between HV(i,j) - HV(i,j+1)
HEF = reshape((1:nrows*(ncols-1)), nrows, ncols-1) + length(VEF(:));

nvars = npixels;
nfac = length(VEF(:)) + length(HEF(:));

G = zeros(nvars, nfac);
N = length(ns);
eclass = zeros(1, nfac); % eclass(i)=j  means factor i gets its params from factors{j}
vfactor_ndx = 1; % all vertcial edges get their params from factors{1}
hfactor_ndx = 2; % all vertcial edges get their params from factors{2}
for i=1:nrows
  for j=1:ncols
    if i < nrows
      G(H(i:i+1,j), VEF(i,j)) = 1;
      eclass(VEF(i,j)) = vfactor_ndx;
    end
    if j < ncols
      G(H(i,j:j+1), HEF(i,j)) = 1;
      eclass(HEF(i,j)) = hfactor_ndx;
    end
  end
end


% "kitten raised in cage" prior - more likely to see continguous vertical lines
vfactor = tabular_kernel([K K], softeye(K, 0.9));
hfactor = tabular_kernel([K K], softeye(K, 0.5));
factors = cell(1,2);
factors{vfactor_ndx} = vfactor;
factors{hfactor_ndx} = hfactor;

ev_eclass = ones(1,N); % every observation factor gets is params from ofactor
ns = K*ones(1,nvars);
%fg = mk_fgraph_given_ev(G, ns, factors, {ofactor}, num2cell(img), 'equiv_class', eclass, 'ev_equiv_class', ev_eclass);
fg = mk_fgraph_given_ev(G, ns, factors, {ofactor}, img, 'equiv_class', eclass, 'ev_equiv_class', ev_eclass);

bnet2 = fgraph_to_bnet(fg);

% inference


maximize = 1;

engine = {};
engine{1} = belprop_fg_inf_engine(fg, 'max_iter', npixels*2);
engine{2} = jtree_inf_engine(bnet2);
nengines = length(engine);

% on fg, we have already included the evidence
evidence = cell(1,npixels);
tic; [engine{1}, ll(1)] = enter_evidence(engine{1}, evidence, 'maximize', maximize); toc


% on bnet2, we must add evidence to the dummy nodes 
V = fg.nvars;
dummy = V+1:V+fg.nfactors;
N = max(dummy);
evidence = cell(1, N);
evidence(dummy) = {1};
tic; [engine{2}, ll(2)] = enter_evidence(engine{2}, evidence); toc


Ihat = zeros(nrows, ncols, nengines);
for e=1:nengines
  for i=1:nrows
    for j=1:ncols
      m = marginal_nodes(engine{e}, H(i,j));
      Ihat(i,j,e) = argmax(m.T)-1;
    end
  end
end
Ihat
