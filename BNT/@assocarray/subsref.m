function val = subsref(A, S)
% SUBSREF Subscript reference for an associative array
% A('foo') will return the value associated with foo.
% If there are multiple identicaly keys, the first match is returned.
% Currently the search is sequential.

i = 1;
while i <= length(A.keys)
  if strcmp(S.subs{1}, A.keys{i})
    val = A.vals{i};
    return;
  end
  i = i + 1;
end
error(['can''t find ' S.subs{1}])
