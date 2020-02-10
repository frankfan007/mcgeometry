function dcsdatastruct=prepare_migraine_dcs_data(dcs_file)
%
% dcsdatastruct=prepare_migraine_dcs_data(dcs_file)
%
% function specifically to process Juliette and Phoebe's migraine data
% input file is ***_SyncData.mat

% author: Melissa Wu, <mwu22@mgh.harvard.edu>
% this function is part of the mcgeometry toolbox,
%(https://github.com/wumelissa/mc_geometry)
%%

load(dcs_file.measurement_file);

%% filtering and creating time and tau arrays

% filter out the bad correlation curves
idx=1;
for tpt=1:size(corrset,1)
    temp_corr=squeeze(corrset(tpt,:,:));
    if max(temp_corr(:))<50 && ~all(temp_corr(:)==0)
        corrset_filter(idx,:,:)=corrset(tpt,:,:);
        corrset_times_filter(idx,:)=corrset_times(tpt,:);
        corrset_intensity_filter(idx,:)=corrset_intensity(tpt,:);
        idx=idx+1;
    end
end

% calculating the elapsed time between each correlation curve
temp_time_arr(1)=0;
for tpt=2:size(corrset_times_filter,1)
    time=etime(corrset_times_filter(tpt,:),corrset_times_filter(tpt-1,:));
    temp_time_arr(tpt)=temp_time_arr(tpt-1)+time;
end

% calculating the taus
first_delay=2e-7; % seconds
for I=1:16
    tau(I) = I*first_delay;
end
for J=1:30
    for I=0:7
        tau(I+(J-1)*8+17) = tau((J-1)*8+16+I)+first_delay*(2^J);
    end
end

%% obtaining time segments from triggers

if isfield(dcs_file,'average_time_segments') && dcs_file.average_time_segments==1
    
    if ~exist('s_allHC','var')
        error('No trigger data available');
    end
    
    trig_indices=find(s_allHC>0);
    nb_trig=length(trig_indices);
    trig_seconds=round(t_nirs(trig_indices));
        
    for segment=1:2:2*length(trig_indices)
        idx=floor(segment/2)+1;
        segment_array(segment,1)=trig_seconds(idx)-60;
        segment_array(segment,2)=trig_seconds(idx);
        segment_array(segment+1,1)=trig_seconds(idx);
        segment_array(segment+1,2)=trig_seconds(idx)+30;
    end
    
    dcs_file.average_time_segments=segment_array;
        
end

%%

dcsdatastruct.g2=corrset_filter;
dcsdatastruct.tau=tau;
dcsdatastruct.counts=corrset_intensity_filter;
dcsdatastruct.t=temp_time_arr;

