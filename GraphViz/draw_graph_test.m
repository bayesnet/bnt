% TEST_LAYOUT       Script to test some bayesian net layouts
%

% Change History :
% Date		Time		Prog	Note
% 13-Apr-2000	10:40 PM	ATC	Created under MATLAB 5.3.1.29215a (R11.1)

% ATC = Ali Taylan Cemgil,
% SNN - University of Nijmegen, Department of Medical Physics and Biophysics
% e-mail : cemgil@mbfys.kun.nl 

%bnet = mk_asia_bnet;
%draw_graph(bnet.dag);

% Make the following network (from Jensen (1996) p84 fig 4.17)
%    1
%  / | \
% 2  3  4
% |  |  |
% 5  6  7
%  \/ \/
%  8   9
% where all arcs point downwards

disp('plot directed')
clf;

N = 9;
dag = zeros(N,N);
dag(1,2)=1; dag(1,3)=1; dag(1,4)=1;
dag(2,5)=1; dag(3,6)=1; dag(4,7)=1;
dag(5,8)=1; dag(6,8)=1; dag(6,9)=1; dag(7,9) = 1;

draw_graph(dag);

pause
clf
disp('plot undirected')
udag = [dag+dag'];
draw_graph(udag);

pause
clf
disp('plot mixed')
mg = [dag];
mg(2,1) = 1; mg(8,5) = 1;
draw_graph(mg);



