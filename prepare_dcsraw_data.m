function dcsdatastruct=prepare_dcsraw_data(measurement_filename,g2freq);
%
% dcsdatastruct=prepare_dcsraw_data(dcs_file);
%
% input:
%   dcs_file: structure with fields
%       measurement_file: name of measurement file to load

% author: Melissa Wu, <mwu22@mgh.harvard.edu>
% this function is part of the mcgeometry toolbox,
%(https://github.com/wumelissa/mc_geometry)
%%

rawDCS=readDCS(measurement_filename);

%%

dcsdatastruct.g2=permute(rawDCS.g2,[3 2 1]);
dcsdatastruct.counts=permute(rawDCS.intensities,[2 1]);
dcsdatastruct.tau=rawDCS.taus';
dcsdatastruct.t=(1/g2freq):(1/g2freq):(size(rawDCS.g2,3)*(1/g2freq));