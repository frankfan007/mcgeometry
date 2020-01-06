function [his_data,photon_indices,photon_fractions_retained,num_tissue_layers]=load_history_file(history_filename,max_photons,varargin)
%
% [his_data,photon_indices,photon_fractions_retained,num_tissue_layers]=load_history_file(history_filename,max_photons,varargin)
%
% reads Monte Carlo photon history file

% input: 
%   history_filename: full filename of photon history file
%   max_photons: 

% output:
%   his_data: photon history array outputted from Monte Carlo simulation
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
%   num_tissue_layers: number of tissue layers of volume, dimension (1,1)
%
% author: Melissa Wu, <mwu22@mgh.harvard.edu>
% this function is part of the mcgeometry toolbox,
%(https://github.com/wumelissa/mc_geometry)
%%

if ~isempty(varargin)
    detector_names=varargin{1};
end

%% reading history file

history_filename=[history_filename '.mch'];
[his_temp,~,~]=loadmch(history_filename);
his_temp=his_temp(:,[1 3:end]);

photon_counts=histcounts(his_temp(:,1),'BinMethod','integer');
photon_counts(photon_counts>max_photons)=max_photons;
photon_indices{1}=[1 photon_counts(1)];

for det_idx=2:length(photon_counts)
    photon_indices{det_idx}=[1+sum(photon_counts(1:(det_idx-1))) sum(photon_counts(1:det_idx))];
end

his_data=zeros(sum(photon_counts),size(his_temp,2));
for det_idx=1:length(photon_counts),
    allphot_idx=find(his_temp(:,1)==(det_idx),max_photons);
    photon_fractions_retained(det_idx)=max_photons/length(find(his_temp(:,1)==det_idx));
    curr_idxs=photon_indices{det_idx}(1):photon_indices{det_idx}(2);
    his_data(curr_idxs,:)=his_temp(allphot_idx,:);
end
    
photon_fractions_retained(photon_fractions_retained>1)=1;
num_tissue_layers=(size(his_data,2)-1)/2;
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
