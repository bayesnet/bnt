function CPD = vonMises_CPD( bnet,self,varargin)
%VONMISES_CPD Makes a conditional linear Von Mises Distribution
%   CPD = vonMises_CPD(bnet, node, ...) will create a CPD with random parameters,
% where node is the number of a node in this equivalence class.

CPD = init_fields;
 
CPD = class(CPD, 'vonMises_CPD', generic_CPD(0));

args = varargin;
ns = bnet.node_sizes;
ps = parents(bnet.dag, self);
dps = myintersect(ps, bnet.dnodes);
cps = myintersect(ps, bnet.cnodes);
fam_sz = ns([ps self]);

CPD.self = self;
CPD.sizes = fam_sz;

% Figure out which (if any) of the parents are discrete, and which cts, and how big they are
% dps = discrete parents, cps = cts parents
CPD.cps = find_equiv_posns(cps, ps); % cts parent index
CPD.dps = find_equiv_posns(dps, ps);
ss = fam_sz(end);
psz = fam_sz(1:end-1);
dpsz = prod(psz(CPD.dps));
cpsz = sum(psz(CPD.cps));

% set default params
CPD.mean = randn(ss, dpsz);
CPD.con = 100*repmat(eye(ss), [1 1 dpsz]);    
CPD.weights = randn(ss, cpsz, dpsz);
CPD.cov_type = 'full';
CPD.tied_cov = 0;
CPD.clamped_mean = 0;
CPD.clamped_cov = 0;
CPD.clamped_weights = 0;
CPD.con_prior_weight = 0.01;
CPD.cov_prior_entropic = 0;
nargs = length(args);
if nargs > 0
  CPD = set_fields(CPD, args{:});
end

% Make sure the matrices have 1 dimension per discrete parent.
% Bug fix due to Xuejing Sun 3/6/01
CPD.mean = myreshape(CPD.mean, [ss ns(dps)]);
CPD.con = myreshape(CPD.con, [ss ss ns(dps)]);
CPD.weights = myreshape(CPD.weights, [ss cpsz ns(dps)]);

% Precompute indices into block structured  matrices
% to speed up CPD_to_lambda_msg and CPD_to_pi
cpsizes = CPD.sizes(CPD.cps);
CPD.cps_block_ndx = cell(1, length(cps));
for i=1:length(cps)
  CPD.cps_block_ndx{i} = block(i, cpsizes);
end

%learning - ESS Statistics
CPD.Wsum = zeros(dpsz,1);
CPD.WYsum = zeros(ss, dpsz);
CPD.WXsum = zeros(cpsz, dpsz);
CPD.WYYsum = zeros(ss, ss, dpsz);
CPD.WXXsum = zeros(cpsz, cpsz, dpsz);
CPD.WXYsum = zeros(cpsz, ss, dpsz);
end

function CPD = init_fields()
% This ensures we define the fields in the same order 
% no matter whether we load an object from a file,
% or create it from scratch. (Matlab requires this.)

CPD.self = [];
CPD.sizes = [];
CPD.cps = [];
CPD.dps = [];
CPD.mean = [];
CPD.con = [];
CPD.weights = [];
CPD.clamped_mean = [];
CPD.clamped_cov = [];
CPD.clamped_weights = [];
CPD.cov_type = [];
CPD.tied_cov = [];
CPD.Wsum = [];
CPD.WYsum = [];
CPD.WXsum = [];
CPD.WYYsum = [];
CPD.WXXsum = [];
CPD.WXYsum = [];
CPD.nsamples = [];
CPD.nparams = [];            
CPD.con_prior_weight = [];
CPD.cov_prior_entropic = [];
CPD.useC = [];
CPD.cps_block_ndx = [];

end
