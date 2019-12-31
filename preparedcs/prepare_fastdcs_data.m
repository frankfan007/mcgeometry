function dcsdatastruct=prepare_fastdcs_data(dcs_file);

% input:
%   dcs_file: structure with fields
%       measurement_file: name of measurement file to load

%%

load(dcs_file.measurement_file);

g2freq=dcs_file.g2freq;
%%

dcsdatastruct.g2=permute(all_g2,[3 2 1]);
dcsdatastruct.counts=all_countrate;
dcsdatastruct.tau=taus;
dcsdatastruct.t=(1/g2freq):(1/g2freq):(size(all_g2,3)*(1/g2freq));
