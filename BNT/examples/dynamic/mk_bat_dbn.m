function [bnet, names] = mk_bat_dbn()
% MK_BAT_DBN Make the BAT DBN
% [bnet, names] = mk_bat_dbn()
% See
% - Forbes, Huang, Kanazawa and Russell, "The BATmobile: Towards a Bayesian Automated Taxi", IJCAI 95
% - Boyen and Koller, "Tractable Inference for Complex Stochastic Processes", UAI98.

names = {'LeftClr', 'RightClr', 'LatAct', 'Xdot', 'InLane', 'FwdAct', ...
      'Ydot', 'Stopped', 'EngStatus', 'FBStatus', ...
      'LeftClrSens', 'RightClrSens', 'TurnSignalSens', 'XdotSens', 'YdotSens', ...
      'FYdotDiffSens', 'FclrSens', 'BXdotSens', 'BclrSens', 'BYdotDiffSens', ...
      'SensorValid', 'FYdotDiff', 'FcloseSlow', 'Fclr', 'BXdot', 'BcloseFast', 'Bclr', 'BYdotDiff'};
ss = length(names);

intrac = {...
      'LeftClr', 'LeftClrSens';
  'RightClr', 'RightClrSens';
  'LatAct', 'TurnSignalSens'; 'LatAct', 'Xdot';
  'Xdot', 'XdotSens';
  'FwdAct', 'Ydot';
  'Ydot', 'YdotSens'; 'Ydot', 'Stopped';
  'EngStatus', 'Ydot'; 'EngStatus', 'FYdotDiff'; 'EngStatus', 'Fclr'; 'EngStatus', 'BXdot';
  'SensorValid', 'XdotSens';   'SensorValid', 'YdotSens';
  'FYdotDiff', 'FYdotDiffSens'; 'FYdotDiff', 'FcloseSlow';
  'FcloseSlow', 'FBStatus';
  'Fclr', 'FclrSens'; 'Fclr', 'FcloseSlow';
  'BXdot', 'BXdotSens';
  'Bclr', 'BclrSens'; 'Bclr', 'BXdot'; 'Bclr', 'BcloseFast';
  'BcloseFast', 'FBStatus';
  'BYdotDiff', 'BYdotDiffSens'; 'BYdotDiff', 'BcloseFast'};
[intra, names] = mk_adj_mat(intrac, names, 1);


interc = {...
      'LeftClr', 'LeftClr'; 'LeftClr', 'LatAct';
  'RightClr', 'RightClr'; 'RightClr', 'LatAct';
  'LatAct', 'LatAct'; 'LatAct', 'FwdAct';
  'Xdot', 'Xdot'; 'Xdot', 'InLane';
  'InLane', 'InLane'; 'InLane', 'LatAct';
  'FwdAct', 'FwdAct';
  'Ydot', 'Ydot';
  'Stopped', 'Stopped';
  'EngStatus', 'EngStatus';
  'FBStatus', 'FwdAct'; 'FBStatus', 'LatAct'};
inter = mk_adj_mat(interc, names, 0);  

obs = {'LeftClrSens', 'RightClrSens', 'TurnSignalSens', 'XdotSens', 'YdotSens', 'FYdotDiffSens', ...
      'FclrSens', 'BXdotSens', 'BclrSens', 'BYdotDiffSens'};

for i=1:length(obs)
  onodes(i) = strmatch(obs{i}, names); %stringmatch(obs{i}, names);
end
onodes = sort(onodes);

dnodes = 1:ss; 
ns = 2*ones(1,ss); % binary nodes
bnet = mk_dbn(intra, inter, ns, 'discrete', dnodes, 'observed', onodes, 'eclass2', (1:ss)+ss);

% make rnd params
for i=1:2*ss
  bnet.CPD{i} = tabular_CPD(bnet, i);
end

