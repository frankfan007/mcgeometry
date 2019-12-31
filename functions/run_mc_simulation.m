function run_mc_simulation(dir_struct,mc_param)

% runs Monte Carlo simulation
% 
% input:
%   dir_struct: structure with fields
%       input_filename: full filename of input file
%       mcx_bin: path containing mcx executable
%   mc_param:
%       gpu_number: number of GPU to use
%       max_detected_photons: maximum number of detected photons

%%

idx=regexp(dir_struct.input_filename,'/','once');
temp_name=dir_struct.input_filename(idx:end);
temp_name=strrep(temp_name,'\','/');
all_words=strsplit(temp_name,'/');

final_filename='';
for word=1:length(all_words)
    final_filename=strcat(final_filename,'/',all_words{word});
end
final_filename=final_filename(2:end);

if ~isunix
    fprintf(['Please run the following command in a Linux workstation:\n'...
        '%s -f %s -m 1 -D P -G %d -H %d\n'],dir_struct.mcx_bin,final_filename,mc_param.gpu_number,mc_param.max_detected_photons);
    uiwait(msgbox('Dismiss this when simulation is complete'));
else
    fprintf(['Running simulation command:\n'...
        '%s -f %s -m 1 -D P -G %d -H %d\n'],dir_struct.mcx_bin,final_filename,mc_param.gpu_number,mc_param.max_detected_photons);   
    system(sprintf('%s -f %s -m 1 -D P -G %d -H %d',dir_struct.mcx_bin,final_filename,mc_param.gpu_number,mc_param.max_detected_photons));
end