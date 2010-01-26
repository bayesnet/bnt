function T = convert_to_table(CPD, domain, evidence)
% CONVERT_TO_TABLE Convert a mlp CPD to a table, incorporating any evidence 
% T = convert_to_table(CPD, domain, evidence)

self = domain(end);                    
ps = domain(1:end-1);                               % self' parents                                       
%cps = myintersect(ps, cnodes);                      % self' continous parents      
cnodes     = domain(CPD.cpndx);
cps        = myintersect(ps, cnodes);
odom = domain(~isemptycell(evidence(domain)));      % obs nodes in the net
assert(myismember(cps, odom));                      % !ALL the CTS parents must be observed!
ns(cps)=1;
dps = mysetdiff(ps, cps);                           % self' discrete parents                                                    
dobs = myintersect(dps, odom);                      % discrete obs parents

% Extract the params compatible with the observations (if any) on the discrete parents (if any)

if ~isempty(dobs),
    dvals = cat(1, evidence{dobs});             
    ns_eff= CPD.sizes;                               % effective node sizes              
    ens=ns_eff;
    ens(dobs) = 1;                              
    S=prod(ens(dps));
    subs = ind2subv(ens(dps), 1:S);
    mask = find_equiv_posns(dobs, dps);        
    for i=1:length(mask),
        subs(:,mask(i)) = dvals(i);
    end     
    support = subv2ind(ns_eff(dps), subs)';
else 
    ns_eff= CPD.sizes;
    support=[1:prod(ns_eff(dps))];
end

W1=[]; b1=[]; W2=[]; b2=[];

W1 = CPD.W1(:,:,support);
b1= CPD.b1(support,:);
W2 = CPD.W2(:,:,support);
b2= CPD.b2(support,:);
ns(odom) = 1;
dpsize = prod(ns(dps));                             % overall size of the self' discrete parents  

x = cat(1, evidence{cps});    
ndata=size(x,2);

if ~isempty(evidence{self})                         %
    app=struct(CPD);                                %
    ns(self)=app.mlp{1}.nout;                       % pump up self to the original dimension if observed
    clear app;                                      %
end                                                 %

T =zeros(dpsize, ns(self));                         %
for i=1:dpsize                                      %                 
    W1app = W1(:,:,i);                              % 
    b1app = b1(i,:);                                % 
    W2app = W2(:,:,i);                              % 
    b2app = b2(i,:);                                % for each of the dpsize combinations of self'parents values 
    z = tanh(x(:)'*W1app + ones(ndata, 1)*b1app);   % we tabulate the corrisponding glm model
    a = z*W2app + ones(ndata, 1)*b2app;             % (element of the cell array CPD.glim)
    appoggio = normalise(exp(a));                   %
    T(i,:)=appoggio;                                %
    W1app=[]; W2app=[]; b1app=[]; b2app=[];         %
    z=[]; a=[]; appoggio=[];                        %
end                                                 %                

if ~isempty(evidence{self})
    appoggio=[];                            %
    appoggio=zeros(1,ns(self));             %
    r = evidence{self};                     %...if self is observed => in output there's only the probability of the 'true' class
    for i=1:dpsize                          % 
          appoggio(i)=T(i,r);               % 
    end
    T=zeros(dpsize,1);
    for i=1:dpsize
        T(i,1)=appoggio(i);                        
    end
    clear appoggio;
    ns(self) = 1;
end
