function dcsdatastruct=prepare_sampledcs_data(dcs_file)
% prepares sample DCS data

load([dcs_file.measurement_file]);

g2freq=dcs_file.g2freq;

%%

dcsdatastruct.g2=g2;
dcsdatastruct.counts=intensities;
dcsdatastruct.tau=taus;
dcsdatastruct.t=(1/g2freq):(1/g2freq):(size(g2,3)*(1/g2freq));