function write_input_file(input_filename,mc_param)

% input:
%   input_filename: full filename of input file
%   mc_param: structure with subfields
%       num_phot_launched: number of photons launched
%       source: array for source position, dimension (1,3)
%       source_unit_vec: vector array for initial direction of photon launch, dimension (1,3)
%       volume_name_unix: filename of volume to be referenced
%       vol: volume that is loaded from mc_param.volume_name
%       num_tissue_layers: number of tissue layers
%       mus: array of scattering values, dimension (1, ntissue layers)
%       mua: single absorption value for all the layers
%       g: anistropy
%       n: refraction index
%       detectors: array of detector positions, dimension (ndetectors, 3)
%       detector_radii: array of detector radii, dimension (1,ndetectors)

%%

% creating file content array
file_content={
    sprintf('%d',mc_param.num_phot_launched);... %
    '4896744';...
    sprintf('%3.3f %3.3f %3.3f', mc_param.source);...
    sprintf('%3.3f %3.3f %3.3f', mc_param.source_unit_vec);...
    '0 9e-09 9e-09';...
    [mc_param.volume_name_unix];...
    sprintf('1 %d 1 %d',size(mc_param.vol,1),size(mc_param.vol,1));...
    sprintf('1 %d 1 %d',size(mc_param.vol,2),size(mc_param.vol,2));...
    sprintf('1 %d 1 %d',size(mc_param.vol,3),size(mc_param.vol,3));...
    sprintf('%d',mc_param.num_tissue_layers)};

idx=length(file_content)+1;

for layer=1:mc_param.num_tissue_layers
    file_content{idx}=sprintf('%1.6f %1.2f %1.2f %1.2f',mc_param.mus(layer)/(1-mc_param.g),mc_param.g,mc_param.mua,mc_param.n);
    idx=idx+1;
end

file_content{idx}=sprintf('%d 1',size(mc_param.detectors,1));
idx=idx+1;

for det=1:size(mc_param.detectors,1)   
    file_content{idx}=sprintf('%3.3f %3.3f %3.3f %1.2f', mc_param.detectors(det,:),mc_param.det_radii(det));
    idx=idx+1;
end

% writing out all the information to the file
fid=fopen(input_filename,'w');

for line=1:length(file_content)
    fprintf(fid,[file_content{line} '\n']);
end

fclose(fid);