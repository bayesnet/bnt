function T = convert_to_table(CPD, domain, evidence)
% CONVERT_TO_TABLE Convert a softmax CPD to a table, incorporating any evidence 
% T = convert_to_table(CPD, domain, evidence)

self       = domain(end);             
ps         = domain(1:end-1);                            
cnodes     = domain(CPD.cpndx);
cps        = myintersect(ps, cnodes);
dps        = domain(CPD.dpndx); 
dps_as_cps = domain(CPD.dps_as_cps.ndx);
all_dps    = union(dps,dps_as_cps);
odom       = domain(~isemptycell(evidence(domain))); 
if ~isempty(cps), assert(myismember(cps, odom)); end % all cts parents must be observed

ns         = zeros(1, max(domain));
ns(domain) = CPD.sizes;
ens        = ns; % effective node sizes
ens(odom)  = 1;

% dpsize >= glimsz because the glm parameters are tied across the dps_as_cps parents
dpsize       = prod(ens(all_dps)); % size of ALL self'discrete parents
dpvals       = cat(1, evidence{myintersect(all_dps, odom)});
cpvals       = cat(1, evidence{cps});
if ~isempty(dps_as_cps),
  separator          = CPD.dps_as_cps.separator;
  dp_as_cpmap        = find_equiv_posns(dps_as_cps, all_dps);
  dops_map           = find_equiv_posns(myintersect(all_dps, odom), all_dps);
  puredp_map         = find_equiv_posns(dps, all_dps);
  subs               = ind2subv(ens(all_dps), 1:prod(ens(all_dps)));
  if ~isempty(dops_map), subs(:,dops_map) = subs(:,dops_map)+repmat(dpvals(:)',[size(subs,1) 1])-1; end
end

[w,b] = extract_params(CPD);
T = zeros(dpsize, ns(self));                                       
for i=1:dpsize,    
  active_glm  = i;
  dp_as_cpvals=zeros(1,sum(ns(dps_as_cps)));                                                                  
  if ~isempty(dps_as_cps),                          
    active_glm = max([1,subv2ind(ns(dps), subs(i,puredp_map))]);
    % Extract the params compatible with the observations (if any) on the 'pure' discrete parents (if any)
    where_one = separator + subs(i,dp_as_cpmap);
    % and get in the dp_as_cp parents...
    dp_as_cpvals(where_one)=1;                    
  end                                               
  T(i,:) = normalise(exp([dp_as_cpvals(:); cpvals(:)]'*w(:,:,active_glm) + b(:,active_glm)'));
end
if myismember(self, odom)
  r = evidence{self};
  T = T(:,r);
end

T = myreshape(T, ens(domain));
