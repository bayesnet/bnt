function stop = is_F2_true_D3(vals)
% function stop = is_F2_true_D3(vals)
% 
% If vals(F2)=2 then level 2 has finished, so we return stop=1
% to stop sample_dbn. Otherwise we return stop=0.
% We assume this is for a D=3 level HHMM.

Q1 = 1; Q2 = 2; Q3 = 3; F3 = 4; F2 = 5; Onode = 6;
stop = 0;
if (iscell(vals) & vals{F2}==2) | (~iscell(vals) & vals(F2)==2)
  stop = 1;
end
