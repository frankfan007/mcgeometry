function rawDCS=readDCS(filename)
%
% rawDCS=readDCS(filename)
%
% reads DCS file from MetaOx
% input:
%   filename: full filename of .dcsraw file

% author: Melissa Wu, <mwu22@mgh.harvard.edu>
% this function is part of the mcgeometry toolbox,
%(https://github.com/wumelissa/mc_geometry)
%%

data=importdata(fullfile(filename));
total_data=data.data;
rowheads=char(data.rowheaders);
tp_indices=find(rowheads(:,1)=='C');

tau_length=tp_indices(2)-tp_indices(1)-1;
num_detects=size(total_data,2);

g2_data=nan(tau_length,num_detects,length(tp_indices));
taus=nan(tau_length,1);
intensities=nan(num_detects,length(tp_indices));

for idx=1:numel(tp_indices)
    elem=tp_indices(idx)+1;
    for det=1:num_detects
        g2_data(:,det,idx)=total_data(elem:elem+tau_length-1,det);
        intensities(det,idx)=total_data(elem-1,det);
    end
    for I=2:(tau_length+1)
        taus(I-1)=str2num(rowheads(I,:));
    end
end

rawDCS.g2=g2_data;
rawDCS.taus=taus;
rawDCS.intensities=intensities;

        