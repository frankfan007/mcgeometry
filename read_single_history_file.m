function [photon_data,photon_indices,photon_fractions_retained]=read_single_history_file(filename,num_layers,max_photons, varargin)
% reads Monte Carlo photon history file from tMCimg

% inputs:
%   filename: name of .his file WITHOUT the extension
%   num_layers: number of layers in binary volume used for MC simulation
%   max_photons: maximum number of photons to save per detector
% output:
%   photon_data: photon history array outputted from Monte Carlo simulation
%       array with path length and momentum transfer information stored for each detected photon
%       dimension (number of detected photons, 1 + 2*(number of tissue layers))
%       column 1 stores the index of the detector that the photon hit
%       columns 2:(1+number of tissue layers) stores the path lengths of the photons for each tissue layer
%       columns (2+number of tissue layers):end stores the momentum transfers of the photons for each tissue layer
%   photon_indices: cell with reference indices for his_array for photons detected by each detector
%       dimension (1, number of detectors)
%           ex: {1:1700,1701:2500,2501:3000}  denotes that columns (photons) 1:1700 in his_array are from detector 1,
%           columns (photons) 1701:2500 in his_array are from detector 2,
%           and columns (photons) 2501:3000 in his_array are from detector 3
%   photon_fractions_retained: fraction of photons saved out of photons detected for each detector
%       dimension (1, number of detectors)

% author: Melissa Wu, <mwu22@mgh.harvard.edu>
% this function is part of the mcgeometry toolbox,
%(https://github.com/wumelissa/mc_geometry)

%%

if ~isempty(varargin)
    detector_names=varargin{1};
end

%%
fid=fopen(filename,'rb');
his_temp=fread(fid,'float');
fclose(fid);

columns=1+num_layers*2;
num_photons=length(his_temp)/columns;
his_temp=transpose(reshape(his_temp,columns,num_photons));
photon_counts=histcounts(his_temp(:,1),'BinMethod','integer');
photon_counts(photon_counts>max_photons)=max_photons;
photon_indices{1}=[1 photon_counts(1)];
for det_idx=2:length(photon_counts)
    photon_indices{det_idx}=[1+sum(photon_counts(1:(det_idx-1))) sum(photon_counts(1:det_idx))];
end

photon_data=zeros(sum(photon_counts),columns);
for det_idx=1:length(photon_counts),
    allphot_idx=find(his_temp(:,1)==(det_idx-1),max_photons);
    photon_fractions_retained(det_idx)=max_photons/length(find(his_temp(:,1)==(det_idx-1)));
    curr_idxs=photon_indices{det_idx}(1):photon_indices{det_idx}(2);
    photon_data(curr_idxs,:)=his_temp(allphot_idx,:);
end

%% print box

if exist('detector_names','var')
    for det=1:length(photon_fractions_retained)
        text_cell{det}=sprintf('%s detector: %2.0f percent of photons retained',detector_names{det},photon_fractions_retained(det)*100);
    end
else
    for det=1:length(photon_fractions_retained)
        text_cell{det}=sprintf('Detector %d: %2.0f percent of photons retained',det,photon_fractions_retained(det)*100);
    end
end

print_box_with_height(text_cell,100);
