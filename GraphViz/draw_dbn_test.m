% TEST_LAYOUT_DBN		Script to test some DBN layouts
%

% Change History :
% Date		Time		Prog	Note
% 17-Apr-2000	 2:40 PM	ATC	Created under MATLAB 5.3.1.29215a (R11.1)

% ATC = Ali Taylan Cemgil,
% SNN - University of Nijmegen, Department of Medical Physics and Biophysics
% e-mail : cemgil@mbfys.kun.nl 

disp('draw mhmm1')
clf
set(gcf, 'pos', [0   0   1024   600]);

intra = zeros(3);
intra(1,[2 3]) = 1;
intra(2,3) = 1;
inter = zeros(3);
inter(1,1) = 1;
n = 3;
dnodes = [1 2];
isbox = zeros(n,1); isbox(dnodes) = 1;
unfold = 4;
draw_dbn(intra, inter, 0, unfold, {'Q', 'M', 'Y'}, isbox);


pause
clf
disp('draw water1')
bnet = mk_water_dbn;
unfold = 3;
flip = 1;
[dummyx, dummyy, h] = draw_dbn(bnet.intra, bnet.inter, flip, unfold);

col = rand(size(h,1),3);
for i=1:length(h),
  col = rand(1,3);
  % patches
  set(h(i,2),'facecolor', col); drawnow;
  % text
  set(h(i,1),'color', 1-col); drawnow;
end;


pause
clf
disp('draw BAT static')
% This requires BNT

%sz = get(0, 'ScreenSize');
%sv = get(gcf, 'pos');
%set(gcf, 'units','pix','pos', sz);
[bnet, names] = mk_bat_dbn;
N = size(bnet.intra,1);
G = [bnet.intra bnet.inter;
     zeros(N,N) bnet.intra];
names2 = names;
for i=1:N
  names2{i} = sprintf('%s(1)', names{i});
  names2{i+N} = sprintf('%s(2)', names{i});
end
draw_graph(G, names2);

pause
disp('draw bat dbn')
clf
draw_dbn(bnet.intra, bnet.inter, 0, 2, names2);
