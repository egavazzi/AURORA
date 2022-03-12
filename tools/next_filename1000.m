function [next_filename] = next_filename1000(filename)
% next_filename1000 - Backup files to prevent overwriting them
%   
% Calling:
%   [status, message] = backup1000filesversions(FILENAME)
% Input:
%   FILENAME - name of file to be moved. The file FILENAME will be
%              moved to FILENAME.N, so that one can save file
%              FILENAME multiple times and get the history of the
%              previous files as FILENAME.1, FILENAME.2, ...
%              FILENAME.N 
% Output:
%  status  -  logical scalar, defining the outcome of COPYFILE.
%             1 : copyfile executed successfully. 0 : an error
%             occurred.
%  message - string, defining the error or warning message.
%            string : copyfile executed successfully. message
%            : an error or warning message, as applicable.

%   Copyright © 2011 M V, Bjorn Gustavsson <bjorn.gustavsson@irf.se>, 
%   This is free software, licensed under GNU GPL version 2 or later

sprintf('\nBacking up the old file as...')
counter = '001';
while exist([filename, '.', counter],'file')
  counter = num2str(sprintf('%03d', str2num(counter)+1));
end
disp(sprintf('%s', [filename, '.', counter]))
[status, message] = movefile(filename, [filename, '.', counter]);
