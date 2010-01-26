function risultati = fhme(net, nodes_info, data, n)
%HMEFWD	Forward propagation through an HME model
%
% Each row of the (n x class_num) matrix 'risultati' containes the estimated class posterior prob.
%
% ----------------------------------------------------------------------------------------------------
% -> pierpaolo_b@hotmail.com   or   -> pampo@interfree.it
% ----------------------------------------------------------------------------------------------------
%
ns=net.node_sizes;
if nargin==3
    ndata=n;
else
    ndata=size(data, 1);
end
altezza=size(ns,2);
coeff=cell(altezza-1,1);
for m=1:ndata
    %- i=2 --------------------------------------------------------------------------------------
    s=struct(net.CPD{2});    
    if nodes_info(1,2)==0,
        mu=[]; W=[]; predict=[];
        mu=s.mean(:,:);
        W=s.weights(:,:,:);
        predict=mu(:,:)+W(:,:,:)*data(m,:)';            
        coeff{1,1}=predict';            
    elseif nodes_info(1,2)==1,
        coeff{1,1}=fglm(s.glim{1}, data(m,:));
    else,
        coeff{1,1}=fmlp(s.mlp{1}, data(m,:));
    end
    %----------------------------------------------------------------------------------------------
    if altezza>3,
        for i=3:altezza-1,
            s=[]; f=[]; dpsz=[];
            f=family(net.dag,i); f=f(2:end-1); dpsz=prod(ns(f));
            s=struct(net.CPD{i});
            for j=1:dpsz,
                if nodes_info(1,i)==1,
                    coeff{i-1,1}(j,:)=coeff{i-2,1}(1,j)*fglm(s.glim{j}, data(m,:));
                else
                    coeff{i-1,1}(j,:)=coeff{i-2,1}(1,j)*fmlp(s.mlp{j}, data(m,:));
                end
            end       
            app=cat(2, coeff{i-1,1}(:)); coeff{i-1,1}=app'; clear app;
        end
    end
    %- i=altezza ----------------------------------------------------------------------------------
    if altezza>2,
        i=altezza;
        s=[]; f=[]; dpsz=[];
        f=family(net.dag,i); f=f(2:end-1); dpsz=prod(ns(f));
        s=struct(net.CPD{i});
        if nodes_info(1,i)==0,            
            mu=[]; W=[];
            mu=s.mean(:,:);
            W=s.weights(:,:,:);
        end
        for j=1:dpsz,
            if nodes_info(1,i)==0,            
                predict=[];
                predict=mu(:,j)+W(:,:,j)*data(m,:)';            
                coeff{i-1,1}(j,:)=coeff{i-2,1}(1,j)*predict';            
            elseif nodes_info(1,i)==1,
                coeff{i-1,1}(j,:)=coeff{i-2,1}(1,j)*fglm(s.glim{j}, data(m,:));
            else
                coeff{i-1,1}(j,:)=coeff{i-2,1}(1,j)*fmlp(s.mlp{j}, data(m,:));
            end
        end
    end
    %----------------------------------------------------------------------------------------------
    risultati(m,:)=sum(coeff{altezza-1,1},1);
    clear coeff; coeff=cell(altezza-1,1);
end
return

%-------------------------------------------------------------------

function [y, a] = fglm(net, x)
%GLMFWD	Forward propagation through 1-layer net->GLM statistical model

ndata = size(x, 1);

a = x*net.w1 + ones(ndata, 1)*net.b1;

nout = size(a,2);
% Ensure that sum(exp(a), 2) does not overflow
maxcut = log(realmax) - log(nout);
% Ensure that exp(a) > 0
mincut = log(realmin);
a = min(a, maxcut);
a = max(a, mincut);
temp = exp(a);
y = temp./(sum(temp, 2)*ones(1,nout));

%-------------------------------------------------------------------

function [y, z, a] = fmlp(net, x)
%MLPFWD	Forward propagation through 2-layer network.

ndata = size(x, 1);

z = tanh(x*net.w1 + ones(ndata, 1)*net.b1);
a = z*net.w2 + ones(ndata, 1)*net.b2;  
temp = exp(a);
nout = size(a,2);
y = temp./(sum(temp,2)*ones(1,nout));

%-------------------------------------------------------------------
