function f = fullfileKPM(varargin)
% fullfileKPM Concatenate strings with file separator, then convert it to a/b/c
% function f = fullfileKPM(varargin)

f = fullfile(varargin{:});
f = strrep(f, '\', '/');
