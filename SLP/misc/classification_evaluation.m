function [ratio, ratiominus, ratioplus, proba_post, yt] = classification_evaluation(bnet, BDT, class)
% Computes the classification ratio of a bnet structure on a test dataset BDT
% [ratio ratiominus ratioplus] = classification_evaluation(bnet, BDT, class)
% 
% [ratiominus rationplus] is the 95 percent confident interval.
% results are in percentage [0 100].
%
%  francois.olivier.c.h@gmail.com
%

    [proba_post,engine] = inference(bnet, mat_to_bnt(BDT), class);
               [tmp yt] = max(proba_post, [],2);
                  [N L] = size(BDT);
                  count = length(find(BDT(class,:)==yt'));
                  ratio = 100*count/L;
[ratiominus, ratioplus] = confiance(ratio,L);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [I, J] = confiance(t, N)
% Compute the 95 percent confident interval
%
% see Y. Bennani and F. Bossaert, 
%     Predictive neural networks for traffic disturbance detection in the telephone network
%     In Proceedings of IMACS-CESA 1996, Lille, France.

Z    = 1.96; % this value for the 95 percent confident interval
T    = t/100;
tmp  = (Z*Z)/N;
D    = 1+tmp;
N1   = T+tmp/2;
tmp2 = T*(1-T)/N + tmp/(4*N);
N2   = Z*sqrt(tmp2);
I    = 100*(N1-N2)/D;
J    = 100*(N1+N2)/D;
