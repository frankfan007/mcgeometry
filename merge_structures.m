function finalstruct=merge_structures(structure1,structure2)
%
% finalstruct=merge_structures(structure1,structure2)
%
% function to merge two structures together

% author: Melissa Wu, <mwu22@mgh.harvard.edu>
% this function is part of the mcgeometry toolbox,
%(https://github.com/wumelissa/mc_geometry)
%%

finalstruct=structure1;
f=fieldnames(structure2);

for I=1:length(f)
    finalstruct.(f{I})=structure2.(f{I});
end