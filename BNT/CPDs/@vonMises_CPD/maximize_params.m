function CPD = maximize_params(CPD, temp)
% MAXIMIZE_PARAMS Set the params of a CPD to their ML values (Von Mises)
% CPD = maximize_params(CPD, temperature)
%
% Temperature is currently ignored.

if ~adjustable_CPD(CPD), return; end


if CPD.clamped_mean
  cl_mean = CPD.mean;
else
  cl_mean = [];
end

if CPD.clamped_cov
  cl_con = CPD.con;
else
  cl_con = [];
end

if CPD.clamped_weights
  cl_weights = CPD.weights;
else
  cl_weights = [];
end

[ssz psz Q] = size(CPD.weights);

[ss cpsz dpsz] = size(CPD.weights); % ss = self size = ssz
if cpsz > CPD.nsamples
  fprintf('gaussian_CPD/maximize_params: warning: input dimension (%d) > nsamples (%d)\n', ...
	  cpsz, CPD.nsamples);
end

prior =  repmat(CPD.con_prior_weight*eye(ssz,ssz), [1 1 Q]);


[CPD.mean, CPD.con, CPD.weights] = ...
    clvm_Mstep(CPD.Wsum, CPD.WYsum, CPD.WYYsum, [], CPD.WXsum, CPD.WXXsum, CPD.WXYsum, ...
	      'cov_type', CPD.cov_type, 'clamped_mean', cl_mean, ...
	      'clamped_cov', cl_con, 'clamped_weights', cl_weights, ...
	      'tied_cov', CPD.tied_cov, ...
	      'cov_prior', prior);

if 0
CPD.mean = reshape(CPD.mean, [ss dpsz]);
CPD.con = reshape(CPD.con, [ss ss dpsz]);
CPD.weights = reshape(CPD.weights, [ss cpsz dpsz]);
end

% Bug fix 11 May 2003 KPM
% clg_Mstep collapses all discrete parents into one mega-node
% but convert_to_CPT needs access to each parent separately
sz = CPD.sizes;
ss = sz(end);

% Bug fix KPM 20 May 2003: 
cpsz = sum(sz(CPD.cps));
%if isempty(CPD.cps)
%  cpsz = 0;
%else
%  cpsz = sz(CPD.cps);
%end
dpsz = sz(CPD.dps);
CPD.mean = myreshape(CPD.mean, [ss dpsz]);
CPD.con = myreshape(CPD.con, [ss ss dpsz]);
CPD.weights = myreshape(CPD.weights, [ss cpsz dpsz]);
