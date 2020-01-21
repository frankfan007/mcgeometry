% DCS file information
dcs_file.dcsraw=0;
dcs_file.fastdcs=1;
dcs_file.filename='concatenated_pressmod_files.mat'; % EDIT
dcs_file.g2freq=1;
dcs_file.det_averaging={1:1,2:4}; 

% average data information
% the following are parameters for averaging the timecourse, and will be processed in the following order:
% decimate first, then moving mean, then averaging
dcs_file.decimate_factor=1; 
dcs_file.moving_mean_window_length=1; 
dcs_file.avg_span=3; 