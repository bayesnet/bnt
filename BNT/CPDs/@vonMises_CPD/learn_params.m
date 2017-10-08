function CPD = learn_params(CPD,fam,data,ns,cnodes)
%LEARN_PARAMS Summary of this function goes here
%   Detailed explanation goes here

%get discrete parent nodes
dpnodes = fam(fam~=CPD.self);
%get parent data
parentdata = cell2mat(data(dpnodes,:));
%get size of parents
class = ns(dpnodes);

num_mixtures = 1;

%calculate numer of mixtures
for i=1:length(class)
    num_mixtures = class(i)*num_mixtures;
end

%calculate combinations for mixtures
combinations = zeros(length(class),num_mixtures);
for i=1:length(class)
   combinations(i,:) = repmat(repelem(1:class(i),2^(i-1)),1,num_mixtures/(class(i)*2^(i-1)));
end

%loop through all the combinations
for mixture=1:length(combinations)
    index = zeros(length(dpnodes),length(parentdata));
    %get current combination
    comb = combinations(:,mixture);
    indices = ones(1,length(parentdata));
    
    %get indices of the evidence for the current combination of parent
    %nodes.
    for i=1:length(dpnodes)
        index(i,:) = parentdata(i,:)==comb(i);
        %uses the and function to check the intesection between the indices
        indices = and(indices,index(i,:));
    end
    %check if there is evidence for this specific mixture
    if any(indices)
        evidence = cell2mat(data(CPD.self,:));
        evidence = evidence(indices);

        n = length(evidence);
        %calculate z
        z=sum(exp(1i*(evidence)))/n;

        %calculate mu
        CPD.mean(mixture) = angle(z);

        %calculate R^2
        Rsqr = ((sum(cos(evidence))*(1/n))^2 + (sum(sin(evidence))*(1/n))^2);

        %calc Re^2
        ReSqr = (n/(n-1))*(Rsqr-(1/n));

        % syms m
        % k = vpasolve(besseli(1,m)==besseli(0,m)*ReSqr,m);
        if mean(evidence)~=evidence(1) && (ReSqr>=0 && ReSqr<1-0.001)
            CPD.con(mixture) = fzero(@(l) (besseli(1,l)/(besseli(0,l))-ReSqr),[0,100]);
        else
            %make it gaussian if not able to solve for k
            CPD.con(mixture) = 100;
        end

    end
end

