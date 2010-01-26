function n = letter2num(l)

% map a-z to 1:26 and A-Z to 27:52
punct_code = [32:47 58:64 91:96 123:126];
digits_code = 48:57;
upper_code = 65:90;
lower_code = 97:122;

c = double(l);
n = c-96;
ndx = find(n <= 0); % upper case
n(ndx) = c(ndx)  - 64 + 26;
