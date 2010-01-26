function montageKPM3(data)
% data{f}(y,x,b) - each frame can have a different size (can can even be empty)

data2 = cell2matPad(data);
montageKPM2(data2)
