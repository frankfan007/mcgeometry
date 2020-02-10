function input_file_status=return_input_file_status(input_filename);
%
% input_file_status=return_input_file_status(input_filename);
%
% returns input file status
%
% input:
%   input_filename: full filename of input file

% author: Melissa Wu, <mwu22@mgh.harvard.edu>
% this function is part of the mcgeometry toolbox,
%(https://github.com/wumelissa/mc_geometry)

%% checking file existence

if isempty(input_filename) || exist(input_filename,'dir')
    error('Please enter string for input_filename\n');
end

if ~exist(input_filename,'file')
    input_file_status=-1;
    subtext='creating input file...';
    print_box('Input file does not exist',subtext,70);
    return
end

if ~exist([input_filename '.mch'],'file')
    input_file_status=0;
    subtext='running simulation...';
    print_box('Input file exists, MC simulation not yet run',subtext,70);
    return
end

input_file_status=1;

%%