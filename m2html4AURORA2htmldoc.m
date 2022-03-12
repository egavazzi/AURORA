%% Script run to generate cross-referenced html-documentation of code.
% 
% This script is run to generate the cross-referenced
% html-documentation of the m-source code in AIDA-tools. This
% requires the m2html package: 
%  http://www.artefact.tk/software/matlab/m2html
% 
% Further, to generate the final version, some additional editing
% of the menu.html files to sort the listings from alphabetical to
% logigal order.

%% Set-up
% The directories to include in the documentation:
dirs4doc = {'.','Data_electron','Data_neutrals','Demo-data',...
            'e_spectras','e_transport_2s', ...
            'e_transport_MS','Flickering_examples',...
            'Plotting','Steady_state_examples','sup-scripts',...
            'tests','tools'};
dirs4doc = {'.','Data_electron','Data_neutrals','Demo-data',...
            'e_spectras','e_transport_2s', ...
            'e_transport_MS',...
            'Plotting','sup-scripts',...
            'tests','tools','utilities'};

%% Back-uping
% I make no automatic back-up attempt here, it might be worthwhile
% to save away the current version - it will, if nothing else, save
% a lot of time when editing the menu.html files.

%% Run
m2html('mfiles',dirs4doc,...
       'htmldir','Html-D20190529',... % Default directory for documentation
       'recursive','on',...         % recurse down into sub-directories
       'global','on',...            % build cross-reference links across AIDA-tools
       'todo','on',...              % Extract TODO for the case that such comments exist
       'graph','on',...             % Make cross-reference call graph
       'template','frame',...       % Style to use
       'index','menu',...           % name for menu-file, goes with frame
       'save','on')                 % Save the global cross-reference calls

%% Make global call-graph:
%
% With m2html a global call-graph can be made with a second
% step. From the m2html package FAQ:
%
% Question: I would like to create a full dependency graph of my
%   functions, m2html creates one for each separate directory. 
% Answer: To build the dependency graph of all your files (by
%   default, it builds a dependency graph for each directory), the
%   solution is: 
%        run m2html recursively using these options:
%        >> m2html(..., 'recursive','on', 'global','on', 'save','on', ...)
%
%        This will parse all the files recursively ('recursive'),
%        will look for all dependencies ('global') and will save
%        the parsing in a MAT-file ('save'). 
%        then you should have a 'm2html.mat' file in the output
%        directory that can be used with 'mdot.m' to create the
%        graph you want:
%
%        >> mdot('PATHTODOC/m2html.mat','m2html.dot');
%        >> !dot -Tpng m2html.dot -o m2html.png 
%
%        and then you have your full dependency graph m2html.png.

