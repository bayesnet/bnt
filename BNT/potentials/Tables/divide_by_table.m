function bigT = divide_by_table(bigT, bigdom, bigsz, smallT, smalldom, smallsz)
% DIVIDE_BY_TABLE 
% bigT = divide_by_table(bigT, bigdom, bigsz, smallT, smalldom, smallsz)
%


Ts = extend_domain_table(smallT, smalldom, smallsz, bigdom, bigsz);
% Replace 0s by 1s before dividing. This is valid, Ts(i)=0 iff Tbig(i)=0.
Ts = Ts + (Ts==0);
%Tbig.T(:) = Tbig.T(:) ./ Ts(:);
bigT(:) = bigT(:) ./ Ts(:);

