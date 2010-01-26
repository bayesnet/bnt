function l = num2letter(n)

% map 1:26 to a-z and 27:52 to A-Z
punct_code = [32:47 58:64 91:96 123:126];
digits_code = 48:57;
upper_code = 65:90;
lower_code = 97:122;

letters = [char(lower_code) char(upper_code)];
l = letters(n);
