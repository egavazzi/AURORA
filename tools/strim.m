function s1 = strim(s)
%STRIM(S) strip the trailing and leading blanks of a string.
%         STRIM(S)is an extension of the MATLAB function DEBLANK(S)
%         and removes both leading and trailing blanks from the string S.
%
% Date 11-08-98
% Written by Stefan Baunack. Please send any bug
% reports or other comments to: s.baunack@ifw-dresden.de.

if ~ischar(s)
	error('Input must be a string.')
end
s1 = deblank(s);
s1 = fliplr(s1);
s1 = deblank(s1);
s1 = fliplr(s1);
