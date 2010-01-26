function CPD = maximize_params(CPD, temp)
% MAXIMIZE_PARAMS Set the params of a CPD to their ML values (dsoftmax) using IRLS
% CPD = maximize_params(CPD, temperature)
% temperature parameter is ignored

% Written by Pierpaolo Brutti

if ~adjustable_CPD(CPD), return; end
options = foptions;

if CPD.verbose
  options(1) = 1;
else
  options(1) = -1;
end
%options(1) = CPD.verbose;

options(2) = CPD.wthresh;
options(3) = CPD.llthresh;
options(5) = CPD.approx_hess;
options(14) = CPD.max_iter;

dpsize = size(CPD.self_vals,3);
for i=1:dpsize,
  mask=find(CPD.eso_weights(:,:,i)>0); % for adapting the parameters we use only positive weighted example
  if  ~isempty(mask),
    if ~isempty(CPD.dps_as_cps.ndx),
        puredp_map = find_equiv_posns(CPD.dpndx, union(CPD.dpndx, CPD.dps_as_cps.ndx)); % find the glm  structure
        subs       = ind2subv(CPD.sizes(union(CPD.dpndx, CPD.dps_as_cps.ndx)),i);       % that corrisponds to the
        active_glm = max([1,subv2ind(CPD.sizes(CPD.dpndx), subs(puredp_map))]);         % i-th 'fictitious' example
        
        CPD.glim{active_glm} = netopt_weighted(CPD.glim{active_glm}, options, CPD.parent_vals(mask',:,i),...
            CPD.self_vals(mask',:,i), CPD.eso_weights(mask',:,i), 'scg');
    else
        alfa = 0.4; if CPD.solo, alfa = 1; end % learning step = 1 <=> self is all alone in the net
        CPD.glim{i} = glmtrain_weighted(CPD.glim{i}, options, CPD.parent_vals(mask',:),...
            CPD.self_vals(mask',:,i), CPD.eso_weights(mask',:,i), alfa);
    end               
  end
  mask=[];
end
