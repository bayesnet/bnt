function data2 = cell2matPad(data)
% data{f}(y,x,b) - each frame can have a different size (can can even be empty)
% data2(y,x,b,f) = zero padded version

Nframes = length(data);
Nbands = -inf;
nr = -inf; nc = -inf;
for f=1:Nframes
  if isempty(data{f}), continue; end
  nr = max(nr, size(data{f},1));
  nc = max(nc, size(data{f},2));
  Nbands = max(Nbands, size(data{f},3));
end    
data2 = zeros(nr, nc, Nbands, Nframes);
for f=1:Nframes
  if isempty(data{f}), continue; end
  data2(1:size(data{f},1), 1:size(data{f},2), :, f) = data{f};
end
if Nbands == 1
  data2 = squeeze(data2); % reshape(data2, [size(data2,1), size(data2,2), Nframes]);
end

