function CPD = update_ess(CPD, fmarginal, evidence, ns, cnodes, hidden_bitv)
% UPDATE_ESS Update the Expected Sufficient Statistics of a CPD (MLP)
% CPD = update_ess(CPD, family_marginal, evidence, node_sizes, cnodes, hidden_bitv)
%
% fmarginal = overall posterior distribution of self and its parents
% fmarginal(i1,i2...,ik,s)=prob(Pa1=i1,...,Pak=ik, self=s| X)
% 
% => 1) prob(self|Pa1,...,Pak)=fmarginal/prob(Pa1,...,Pak) with prob(Pa1,...,Pak)=sum{s,fmarginal}
%       [self estimation -> CPD.self_vals]
% 	  2) prob(Pa1,...,Pak) [SCG weights -> CPD.eso_weights]
%
% Hidden_bitv is ignored

% Written by Pierpaolo Brutti

if ~adjustable_CPD(CPD), return; end

dom = fmarginal.domain;                              
cdom = myintersect(dom, cnodes);                     
assert(~any(isemptycell(evidence(cdom))));           
ns(cdom)=1;

self = dom(end);                                  
ps=dom(1:end-1);                                     
dpdom=mysetdiff(ps,cdom);                            

dnodes = mysetdiff(1:length(ns), cnodes);            

ddom = myintersect(ps, dnodes);                      %
if isempty(evidence{self}),                          % if self is hidden in what follow we must 
    ddom = myintersect(dom, dnodes);                 % consider its dimension
end                                                  % 

odom = dom(~isemptycell(evidence(dom)));    
hdom = dom(isemptycell(evidence(dom)));              % hidden parents in domain
 
dobs = myintersect(ddom, odom);             
dvals = cat(1, evidence{dobs});             
ens = ns;                                            % effective node sizes              
ens(dobs) = 1;                              
                                            
dpsz=prod(ns(dpdom));
S=prod(ens(ddom));
subs = ind2subv(ens(ddom), 1:S);
mask = find_equiv_posns(dobs, ddom);
for i=1:length(mask),
    subs(:,mask(i)) = dvals(i);
end
supportedQs = subv2ind(ns(ddom), subs);

Qarity = prod(ns(ddom));
if isempty(ddom),                      
  Qarity = 1;                         
end                                
fullm.T = zeros(Qarity, 1);
fullm.T(supportedQs) = fmarginal.T(:);

% For dynamic (recurrent) net-------------------------------------------------------------
% ----------------------------------------------------------------------------------------
high=size(evidence,1);                                  % slice height
ss_ns=ns(1:high);                                       % single slice nodes sizes
pos=self;                                               %
slice_num=0;                                            %
while pos>high,                                         % 
    slice_num=slice_num+1;                              % find active slice
    pos=pos-high;                                       % pos=self posistion into a single slice
end                                                     %

last_dim=pos-1;                                         % 
if isempty(evidence{self}),                             % 
    last_dim=pos;                                       %
end                                                     % last_dim=last reshaping dimension      
reg=dom-slice_num*high;
dex=myintersect(reg(find(reg>=0)), [1:last_dim]);       %           
rs_dim=ss_ns(dex);                                      % reshaping dimensions

if slice_num>0,
    act_slice=[]; past_ancest=[];                       %
    act_slice=slice_num*high+[1:high];                  % recover the active slice nodes
    % past_ancest=mysetdiff(ddom, act_slice);
    past_ancest=mysetdiff(ps, act_slice);               % recover ancestors contained into past slices
    app=ns(past_ancest);
    rs_dim=[app(:)' rs_dim(:)'];                        %
end                                                     %
if length(rs_dim)==1, rs_dim=[1 rs_dim]; end            %
if size(rs_dim,1)~=1, rs_dim=rs_dim';    end            %

fullm.T=reshape(fullm.T, rs_dim);                       % reshaping the marginal

% ----------------------------------------------------------------------------------------
% ----------------------------------------------------------------------------------------

% X = cts parent, R = discrete self

% 1) observations vector -> CPD.parents_vals -------------------------------------------------
x = cat(1, evidence{cdom});

% 2) weights vector -> CPD.eso_weights -------------------------------------------------------
if isempty(evidence{self}) % R is hidden
    sum_over=length(rs_dim);
    app=sum(fullm.T, sum_over);    
    pesi=reshape(app,[dpsz,1]);
    clear app;
else
    pesi=reshape(fullm.T,[dpsz,1]);
end

assert(approxeq(sum(pesi),1));

% 3) estimate (if R is hidden) or recover (if R is obs) self'value----------------------------
if isempty(evidence{self})              % R is hidden    
    app=mk_stochastic(fullm.T);         % P(self|Pa1,...,Pak)=fmarginal/prob(Pa1,...,Pak)
    app=reshape(app,[dpsz ns(self)]);   % matrix size: prod{j,ns(Paj)} x ns(self)      
    r=app;
    clear app;
else
    r = zeros(dpsz,ns(self));
    for i=1:dpsz
        if pesi(i)~=0, r(i,evidence{self}) = 1; end
    end
end
for i=1:dpsz
    if pesi(i) ~=0, assert(approxeq(sum(r(i,:)),1)); end
end

CPD.nsamples = CPD.nsamples + 1;            
CPD.parent_vals(CPD.nsamples,:) = x(:)';
for i=1:dpsz
    CPD.eso_weights(CPD.nsamples,:,i)=pesi(i);
    CPD.self_vals(CPD.nsamples,:,i) = r(i,:); 
end
