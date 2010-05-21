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
