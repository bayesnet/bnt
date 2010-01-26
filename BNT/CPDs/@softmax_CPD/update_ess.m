function CPD = update_ess(CPD, fmarginal, evidence, ns, cnodes, hidden_bitv)
% UPDATE_ESS Update the Expected Sufficient Statistics of a softmax node
% function CPD = update_ess(CPD, fmarginal, evidence, ns, cnodes, hidden_bitv)
%
% fmarginal = overall posterior distribution of self and its parents
% fmarginal(i1,i2...,ik,s)=prob(Pa1=i1,...,Pak=ik, self=s| X)
% 
% => 1) prob(self|Pa1,...,Pak)=fmarginal/prob(Pa1,...,Pak) with prob(Pa1,...,Pak)=sum{s,fmarginal}
%       [self estimation -> CPD.self_vals]
% 	  2) prob(Pa1,...,Pak) [WIRLS weights -> CPD.eso_weights]
%
% Hidden_bitv is ignored

% Written by Pierpaolo Brutti

if ~adjustable_CPD(CPD), return; end

domain     = fmarginal.domain;                              
self       = domain(end);          
ps         = domain(1:end-1);                                     
cnodes     = domain(CPD.cpndx);
cps        = myintersect(domain, cnodes);                     
dps        = mysetdiff(ps, cps);                            
dn_use     = dps;
if isempty(evidence{self}) dn_use = [dn_use self]; end % if self is hidden we must consider its dimension  
dps_as_cps = domain(CPD.dps_as_cps.ndx);
odom       = domain(~isemptycell(evidence(domain))); 

ns = zeros(1, max(domain));
ns(domain) = CPD.sizes;     % CPD.sizes = bnet.node_sizes([ps self]);
ens = ns;                   % effective node sizes
ens(odom) = 1;              
dpsize = prod(ns(dps));

% Extract the params compatible with the observations (if any) on the discrete parents (if any)
dops = myintersect(dps, odom);
dpvals = cat(1, evidence{dops});

subs = ind2subv(ens(dn_use), 1:prod(ens(dn_use)));
dpmap = find_equiv_posns(dops, dn_use);
if ~isempty(dpmap), subs(:,dpmap) = subs(:,dpmap)+repmat(dpvals(:)',[size(subs,1) 1])-1; end
supportedQs = subv2ind(ns(dn_use), subs); subs=subs(1:prod(ens(dps)),1:length(dps));
Qarity = prod(ns(dn_use));
if isempty(dn_use), Qarity = 1; end   

fullm.T              = zeros(Qarity, 1);
fullm.T(supportedQs) = fmarginal.T(:);
rs_dim = CPD.sizes;    rs_dim(CPD.cpndx) = 1;           %
if ~isempty(evidence{self}), rs_dim(end)=1; end         % reshaping the marginal
fullm.T              = reshape(fullm.T, rs_dim);        %

% --------------------------------------------------------------------------------UPDATE--

CPD.nsamples = CPD.nsamples + 1;

% 1) observations vector -> CPD.parents_vals ---------------------------------------------
cpvals = cat(1, evidence{cps});

if ~isempty(dps_as_cps),   % ...get in the dp_as_cp parents... 
    separator          = CPD.dps_as_cps.separator;
    dp_as_cpmap        = find_equiv_posns(dps_as_cps, dps);       
    for i=1:dpsize,
        dp_as_cpvals=zeros(1,sum(ns(dps_as_cps)));
        possible_vals = ind2subv(ns(dps),i);
        ll=find(ismember(subs(:,dp_as_cpmap), possible_vals(dp_as_cpmap), 'rows')==1);   
        if ~isempty(ll),
            where_one = separator + possible_vals(dp_as_cpmap);
            dp_as_cpvals(where_one)=1;                            
        end
        CPD.parent_vals(CPD.nsamples,:,i) = [dp_as_cpvals(:); cpvals(:)]';
    end
else
    CPD.parent_vals(CPD.nsamples,:) = cpvals(:)';
end

% 2) weights vector -> CPD.eso_weights ----------------------------------------------------
if isempty(evidence{self}),             % self is hidden
    pesi=reshape(sum(fullm.T, length(rs_dim)),[dpsize,1]);
else
    pesi=reshape(fullm.T,[dpsize,1]);
end
assert(approxeq(sum(pesi),1));          % check

% 3) estimate (if R is hidden) or recover (if R is obs) self'value-------------------------
if isempty(evidence{self})                                  % P(self|Pa1,...,Pak)=fmarginal/prob(Pa1,...,Pak)
    r=reshape(mk_stochastic(fullm.T), [dpsize ns(self)]);   % matrix size: prod{j,ns(Paj)} x ns(self)      
else
    r = zeros(dpsize,ns(self));
    for i=1:dpsize, if pesi(i)~=0, r(i,evidence{self}) = 1; end; end
end
for i=1:dpsize, if pesi(i)~=0, assert(approxeq(sum(r(i,:)),1)); end; end     % check

% 4) save the previous values --------------------------------------------------------------
for i=1:dpsize
    CPD.eso_weights(CPD.nsamples,:,i)=pesi(i);
    CPD.self_vals(CPD.nsamples,:,i) = r(i,:); 
end
