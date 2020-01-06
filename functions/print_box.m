function print_box(text1,subtext,box_width);
%
% print_box(text1,subtext,box_width);
%
% prints a nice text box
% input:
%   text1: title text
%   subtext: subtext
%   box_width: width of box
%
% author: Melissa Wu, <mwu22@mgh.harvard.edu>
% this function is part of the mcgeometry toolbox,
%(https://github.com/wumelissa/mc_geometry)
%%

text_length=length(text1);
stext_length=length(subtext);

fprintf([repmat('=',1,box_width) '\n']);
% fprintf(['|' repmat(' ',1,box_width-2) '|\n']);
fprintf(['|' repmat(' ',1,ceil((box_width-text_length-2)/2)) text1 repmat(' ',1,floor((box_width-text_length-2)/2)) '|\n']);
fprintf(['|' repmat(' ',1,ceil((box_width-stext_length-2)/2)) subtext repmat(' ',1,floor((box_width-stext_length-2)/2)) '|\n']);
% fprintf(['|' repmat(' ',1,box_width-2) '|\n']);
fprintf([repmat('=',1,box_width) '\n']);

