function disp_map_hhmm(bnet)

eclass = bnet.equiv_class;
U = 1; A = 2; C = 3; F = 4;

S = struct(bnet.CPD{eclass(A,2)});
disp('abstract trans')
dispcpt(S.transprob)

S = struct(bnet.CPD{eclass(C,2)});
disp('concrete trans for go left') % UAC AC
dispcpt(squeeze(S.transprob(1,:,:,:,:)))

