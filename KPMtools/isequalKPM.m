function p = isequalKPM(a,b)

if isempty(a) & isempty(b)
  p = 1;
else
  p = isequal(a,b);
end
