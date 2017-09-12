function CPD = learn_params(CPD,fam,data,ns,cnodes)
%LEARN_PARAMS Summary of this function goes here
%   Detailed explanation goes here

%get discrete parent nodes
dpnodes = fam(fam~=CPD.self);
count = 1;
for dpnode=1:length(dpnodes)
    %loop through classes in each parent node
    %get parent data
    parentdata = cell2mat(data(dpnodes(dpnode),:));
    for class=1:ns(dpnodes(dpnode))
       
        %get parent data for a specific class
        index = parentdata==class;
        %get evidence for self
        evidence = cell2mat(data(CPD.self,:));
        %get indexed evidence for self
        evidence = evidence(index);
        n = length(evidence);
        %calculate z
        z=sum(exp(1i*(evidence)))/n;

        %calculate mu
        CPD.mean(count) = angle(z);

        %calculate R^2
        Rsqr = ((sum(cos(evidence))*(1/n))^2 + (sum(sin(evidence))*(1/n))^2);

        %calc Re^2
        ReSqr = (n/(n-1))*(Rsqr-(1/n));

        % syms m
        % k = vpasolve(besseli(1,m)==besseli(0,m)*ReSqr,m);
        if mean(evidence)==evidence(1) || std(evidence)<=0.0708
            CPD.con(count) = 1;
        else
            CPD.con(count) = fzero(@(l) besseli(1,l)-(besseli(0,l)*ReSqr),[0,100]);
        end
        count = count+1;
    end
end

end

