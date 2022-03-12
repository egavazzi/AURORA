function y = in_def2(question,default_value)
% IN_DEF2 - input, with default value
% ===========================================================
% function y = in_def2(question,default_value)
%
% Input		question	Text displayed when asking for input.
%		default_value	This value will be retrurned if
%				the user enters nothing.
%
% Output	The value entered by the user or the default value.
%
% Magnus Gustafsson & Maria Ward 1993
%
% Improved by Peter Rydesaeter 1998
% 
% Used with permission
%============================================================

d = num2str(default_value);

if ( ischar(default_value) )
  
  t = input([question,' {',d,'}: '],'s');
  
else
  
  t = input([question,' {',d(1:length(d)),'}: ']);
  
end

if ( isempty(t) )
  
  y=default_value;
  
else
  
  y=t;
  
end
