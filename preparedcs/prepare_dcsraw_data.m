function dcsdatastruct=prepare_dcsraw_data(dcs_file);

% input:
%   dcs_file: structure with fields
%       measurement_file: name of measurement file to load

%%

rawDCS=readDCS(dcs_file.measurement_file);

g2freq=dcs_file.g2freq;
%%

dcsdatastruct.g2=permute(rawDCS.g2,[3 2 1]);
dcsdatastruct.counts=permute(rawDCS.intensities,[2 1]);
dcsdatastruct.tau=rawDCS.taus';
dcsdatastruct.t=(1/g2freq):(1/g2freq):(size(rawDCS.g2,3)*(1/g2freq));