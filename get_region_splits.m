function regions_to_split=get_region_splits(sup_thicknesses,mid_thicknesses,total_layers)
% creates cell array containing arrays to split layers by
%
% input: 
%   sup_thicknesses: array with various superficial thickness layers
%   mid_thicknesses: array with various middle thickness layers
%   total_layers: total number of tissue layers in volume
%
% output:
%   regions_to_split: cell array with all combinations of superficial and middle thicknesses
%       dimension (length(sup_thicknesses) * length(mid_thicknesses),3)

% author: Melissa Wu, <mwu22@mgh.harvard.edu>
% this function is part of the mcgeometry toolbox,
%(https://github.com/wumelissa/mc_geometry)
%%

if isempty(mid_thicknesses)
    regions_to_split=cell(length(sup_thicknesses),2);
    for t=1:length(sup_thicknesses)
        regions_to_split{t,1}=1:sup_thicknesses(t);
        regions_to_split{t,2}=(sup_thicknesses(t)+1):total_layers;
    end
else
    regions_to_split=cell(length(sup_thicknesses)*length(mid_thicknesses),3);
    for mid=0:(length(mid_thicknesses)-1)
        for t=1:length(sup_thicknesses)
            if sup_thicknesses(t)+mid_thicknesses(mid+1)>=total_layers
                regions_to_split(length(sup_thicknesses)*mid+t,:)={[],[],[]};
            else
                regions_to_split{length(sup_thicknesses)*mid+t,1}=1:(sup_thicknesses(t));
                regions_to_split{length(sup_thicknesses)*mid+t,2}=(sup_thicknesses(t)+1):(sup_thicknesses(t)+mid_thicknesses(mid+1));
                regions_to_split{length(sup_thicknesses)*mid+t,3}=(sup_thicknesses(t)+1+mid_thicknesses(mid+1)):total_layers;
            end
        end
    end
end