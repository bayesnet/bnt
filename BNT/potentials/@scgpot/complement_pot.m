function [margpot, comppot] = complement_pot(pot, keep)
% COMPLEMENT_POT complement means decompose of a potential into its strong marginal and 
% its complement corresponds exactly to the decomposition of a probability distribution 
% into its marginal and conditional
% [margpot, comppot] = complement_pot(pot, keep)

% keep can only include continuous head nodes and discrete nodes
% margpot is the stable CG potential of keep nodes
% comppot is the stable CG potential of others in corresponds exactly to 
% the discomposition of a probability distribution of its marginal and conditional

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of the marginal requires integration over      %
% all variables in csumover. Thus cheadkeep contains all     %
% continuous variables in the marginal potential             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%keyboard;
csumover = mysetdiff(pot.cheaddom, keep);
cheadkeep = mysetdiff(pot.cheaddom, csumover);

nodesizes = zeros(1, max(pot.domain));
nodesizes(pot.ddom) = pot.dsizes;
nodesizes(pot.cheaddom) = pot.cheadsizes;
nodesizes(pot.ctaildom) = pot.ctailsizes;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description of the variables in the marginal domain        %
% For the calculation of a strong marginal first integration %
% over all continuous variables in the head takes place.     %
% The calculation of the marginal over the head variables    %
% might result in a smaller or empty tail                    %
% If there are no head variables, and therefore no tail      %
% variables, left marginalisation over discrete variables    %
% may take place                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
margdom      = mysetdiff(pot.domain,keep);
% margddom   = pot.ddom;
margcheaddom = cheadkeep;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Marginalisation over discrete variables is only allowed when %
% the tail is empty                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
margddom = myintersect(pot.ddom,keep);               % Discrete domain of marginal
margctaildom = myintersect(pot.ctaildom,keep);       % Tail domain
assert(isempty(mysetdiff(pot.ddom,margddom)) | isempty(margctaildom))  


%margctaildom = pot.ctaildom;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Even if marginalisation over continuous variables is only defined %
% for head variables, the marginalisation over haed-variables might %
% result in a smaller tail                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
margctaildom = myintersect(pot.ctaildom,keep);

margcheadsizes = nodesizes(margcheaddom);
margcheadsize = sum(margcheadsizes);
margctailsizes = nodesizes(margctaildom);
margctailsize = sum(margctailsizes);

compdom = pot.domain;
compddom = pot.ddom;
compcheaddom = csumover;
compctaildom = myunion(pot.ctaildom, cheadkeep);
compcheadsizes = nodesizes(compcheaddom);
compcheadsize = sum(compcheadsizes);
compctailsizes = nodesizes(compctaildom);
compctailsize = sum(compctailsizes);

dkeep = myintersect(pot.ddom, keep);
%if dom is only contain discrete node
if isempty(pot.cheaddom)
    dsumover = mysetdiff(pot.ddom, dkeep);
    
    if isempty(dsumover)
        margpot = pot;
        comppot = scgpot([], [], [], []);
        return;
    end
        
    
    I = prod(nodesizes(dkeep));
    J = prod(nodesizes(dsumover));
    sum_map = find_equiv_posns(dsumover, pot.ddom);
    keep_map = find_equiv_posns(dkeep, pot.ddom);
    iv = zeros(1, length(pot.ddom)); % index vector
    p1 = zeros(I,J);
    for i=1:I
        keep_iv = ind2subv(nodesizes(dkeep), i);
        iv(keep_map) = keep_iv;
        for j=1:J
            sum_iv = ind2subv(nodesizes(dsumover), j);
            iv(sum_map) = sum_iv;
            k = subv2ind(nodesizes(pot.ddom), iv);
            potc = struct(pot.scgpotc{k}); % violate object privacy
            p1(i,j) = potc.p;
        end
    end
    p2 = sum(p1,2);
    p2 = p2 + (p2==0)*eps;
    
    margscpot = cell(1, I);
    compscpot = cell(1, I*J);
    iv = zeros(1, length(pot.ddom)); % index vector
    for i=1:I
        margscpot{i} = scgcpot(0, 0, p2(i));
        keep_iv = ind2subv(nodesizes(dkeep), i);
        iv(keep_map) = keep_iv;
        for j=1:J
            sum_iv = ind2subv(nodesizes(dsumover), j);
            iv(sum_map) = sum_iv;
            k = subv2ind(nodesizes(pot.ddom), iv);
            q = p1(i,j)/p2(i);
            compscpot{k} = scgcpot(0, 0, q);
        end
    end
    
    margpot = scgpot(dkeep, [], [], nodesizes, margscpot);
    comppot = scgpot(pot.ddom, [], [], nodesizes,compscpot);
    return;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% head of the potential is not empty %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dsize = pot.dsize;
compscpot = cell(1, dsize);

fmaskh = find_equiv_posns(margcheaddom, compctaildom);
fmaskt = find_equiv_posns(margctaildom, compctaildom);

fh = block(fmaskh, compctailsizes);
ft = block(fmaskt, compctailsizes);


if ~isempty(margcheaddom)
    for i=1:dsize
        potc = struct(pot.scgpotc{i});
        q = 1;
        p = potc.p;
        [A1, A2, B1, B2, C11, C12, C21, C22] = partition_matrix_vec_3(potc.A, potc.B, potc.C, margcheaddom, compcheaddom, nodesizes);

        if ~isempty(margcheaddom)
            margscpot{i} = scgcpot(margcheadsize, margctailsize, p, A1, B1, C11);
        else
            margscpot{i} = scgcpot(margcheadsize, margctailsize, p);
        end 
    
        if ~isempty(compcheaddom)
            if ~isempty(margcheaddom)
                E = A2 - C21*pinv(C11)*A1;
                tmp1 = C21*pinv(C11);
                tmp2 = B2 - C21*pinv(C11)*B1;
                F = zeros(compcheadsize, compctailsize);
                F(:, fh) = tmp1;
                F(:, ft) = tmp2;
                G = C22 - C21*pinv(C11)*C12;
            else
                E = A2;
                F = B2;
                G = C22;
            end
            compscpot{i} = scgcpot(compcheadsize, compctailsize, q, E, F, G);
        else
            compscpot{i} = scgcpot(compcheadsize, 0, q);
        end
        if isempty(margcheaddom)
            margpot = scgpot(margddom, [], [], nodesizes, margscpot);
        else
            margpot = scgpot(margddom, margcheaddom, margctaildom, nodesizes, margscpot);
        end
    end
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Marginalisation took place over all head variables.                               %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate the strong marginal %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    margpot = marginalize_pot(pot,keep);
    mPot    = struct(margpot); 
    for i =1:dsize
        potc = struct(pot.scgpotc{i});  
        % Get the probability of the original potential % 
	q = potc.p;
         
        % Get the configuration defined by the index i%
        config = ind2subv(pot.dsizes,i);
        
        % Calculate the corresponding configuration in the marginal potential
        if isempty(margpot.dsizes)
            % keep == []
	    indMargPot = 1;
        else
            equivPos   = find_equiv_posns(dkeep,pot.ddom);
            indMargPot = subv2ind(margpot.dsizes,config(equivPos));
        end
        % Figure out the corresponding marginal potential
        mPotC = struct(mPot.scgpotc{indMargPot});
        p = mPotC.p;
        if p == 0
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % The following assignment is correct as p is only zero if q is also zero %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            compscpot{i} = scgcpot(compcheadsize,compctailsize,0,potc.A,potc.B,potc.C);
        else
            compscpot{i} = scgcpot(compcheadsize,compctailsize,q/p,potc.A,potc.B,potc.C);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Put all components in one potential %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(compcheaddom)
    comppot = scgpot(compddom, [], [], nodesizes,compscpot);
else
    comppot = scgpot(compddom, compcheaddom, compctaildom, nodesizes,compscpot);
end


