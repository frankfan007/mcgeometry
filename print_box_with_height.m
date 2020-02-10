function print_box_with_height(text_cell,box_width);
%
% print_box_with_height(text_cell,box_width);
%
% prints box with variable lines of text
% input:
%   text_cell: cell with text for each line, dimension (1,nlines)
%   box_width: width of box

% author: Melissa Wu, <mwu22@mgh.harvard.edu>
% this function is part of the mcgeometry toolbox,
%(https://github.com/wumelissa/mc_geometry)
%%

fprintf([repmat('=',1,box_width) '\n']);
for text_idx=1:length(text_cell)
    text1=text_cell{text_idx};
    text_length=length(text1);
    fprintf(['|' repmat(' ',1,ceil((box_width-text_length-2)/2)) text1 repmat(' ',1,floor((box_width-text_length-2)/2)) '|\n']);
end
fprintf([repmat('=',1,box_width) '\n']);