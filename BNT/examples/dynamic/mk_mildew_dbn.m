function bnet = mk_mildew_dbn()

% DBN for foreacasting the gross yield of wheat based on climatic data,
% observations of leaf area index (LAI) and extension of mildew,
% and knowledge of amount of fungicides used and time of usage.
% From Kjaerulff '95.

Fungi=1; Mildew=2; LAI=3; Precip=4; Temp=5; Micro=6; Solar=7; Photo=8; Dry=9;
n = 9;
intra = zeros(n,n);
intra(Mildew, LAI)=1;
intra(LAI,[Micro Photo])=1;
intra(Precip,Micro)=1;
intra(Temp,[Micro Photo])=1;
intra(Solar,Photo)=1;
intra(Photo,Dry)=1;

inter = zeros(n,n);
inter(Fungi,Mildew)=1;
inter(Mildew,Mildew)=1;
inter(LAI,LAI)=1;
inter(Micro,Mildew)=1;
inter(Dry,Dry)=1;

ns = 2*ones(1,n);
bnet = mk_dbn(intra, inter, ns, 'observed', [Photo]);

for e=1:max(bnet.equiv_class(:))
  i = bnet.rep_of_eclass(e);
  bnet.CPD{e} = tabular_CPD(bnet,i);
end
