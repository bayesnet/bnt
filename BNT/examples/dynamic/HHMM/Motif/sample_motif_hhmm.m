%bnet = mk_motif_hhmm('motif_pattern', 'acca', 'background', 't');
bnet = mk_motif_hhmm('motif_pattern', 'accaggggga', 'background', []);

chars = ['a', 'c', 'g', 't'];
Tmax = 100;

for seqi=1:5
  evidence = cell2num(sample_dbn(bnet, 'length', Tmax));
  chars(evidence(end,:))
end
