function MChis=concatenate_layers(his_data,tiss_arr)

% returns concatenated history file with given concatenation parameters and input history file

% input:
%   his_data: photon history array outputted from Monte Carlo simulation
%       array with path length and momentum transfer information stored for each detected photon
%       dimension (number of detected photons, 1 + 2*(number of tissue layers))
%       column 1 stores the index of the detector that the photon hit
%       columns 2:(1+number of tissue layers) stores the path lengths of the photons for each tissue layer
%       columns (2+number of tissue layers):end stores the momentum transfers of the photons for each tissue layer
%   concatenate_tissue_layers: cell array denoting tissue layers to concatenate, dimension (1,number of concatetated layers)
%       example: {1:4,5:13,14:21} would concatenate a 21-layer photon history array into 3 layers;
%       it would concatenate layers 1-4 together, 5-13 together, and 14-21 together

% output:
%   MChis: photon history array that has been concatenated from his_data array
%       dimension (number of detected photons, 1 + 2*(number of concatenated tissue layers))
%       column information would be analogous to information in his_data

%% extracting number of tissue layers

last_cell=tiss_arr{end};
num_tissue_layers=last_cell(end);

%% concatenating tissue layers

pl=nan(size(his_data,1),length(tiss_arr));
mt=nan(size(his_data,1),length(tiss_arr));
for tiss_layer=1:length(tiss_arr)
    pl(:,tiss_layer)=sum(his_data(:,tiss_arr{tiss_layer}+1),2);
    mt(:,tiss_layer)=sum(his_data(:,tiss_arr{tiss_layer}+1+num_tissue_layers),2);
end

MChis=[his_data(:,1) pl mt];
