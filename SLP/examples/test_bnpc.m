close all
clear all

N=5;
names={'A','B','C','D','E'};
dag=zeros(N);
dag(1,2)=1;
dag(2,[3 4])=1;
dag([3 4],5)=1;
[xx yy] = draw_graph(dag,names,ones(N,1));
title('Cheng example.');

node_sizes=2*ones(1,N);
bnet=mk_bnet(dag,node_sizes);

if 1
  disp('Generating DataSet');
  bnet.CPD{1} = tabular_CPD(bnet, 1, [0.4 0.6]);
  bnet.CPD{2} = tabular_CPD(bnet, 2, [0.2 0.3 0.8 0.7]);
  bnet.CPD{3} = tabular_CPD(bnet, 3, [0.1 0.2 0.9 0.8]);
  bnet.CPD{4} = tabular_CPD(bnet, 4, [0.6 0.8 0.4 0.2]);
  bnet.CPD{5} = tabular_CPD(bnet, 5, [0.9 0.8 0.7 0.6 0.1 0.2 0.3 0.4]);
  m=10000;
  cheng = cell(N,m);
  for i=1:m
    cheng(:,i)=sample_bnet(bnet);
  end
  cheng = cell2num(cheng);
  %save -ascii cheng cheng
else
  load -ascii cheng.mat
end

%  profile clear
%  profile on
  [Phase_3, Phase_2, Phase_1, UPhase_3] = learn_struct_bnpc(cheng, node_sizes, 0.05, 0)
%  profile off
%  profile report report_cheng

figure
draw_graph(Phase_1,names,ones(N,1),xx,yy);
title('PhaseI');
figure
draw_graph(Phase_2,names,ones(N,1),xx,yy);
title('PhaseII');
%figure
%draw_graph(UPhase_3,names,ones(N,1),xx,yy);
%title('undirected PhaseIII');
figure
draw_graph(Phase_3,names,ones(N,1),xx,yy);
title('PhaseIII');
