function [g2_data,time_arr,intensities,tau]=pre_process_dcs_data(dcs_file,dcsdatastruct)
%
% [g2_data,time_arr,intensities,tau]=pre_process_dcs_data(dcs_file,dcsdatastruct)
%
% pre-processes DCS data to average detectors and time segments

% input:
%   dcs_file: structure with fields:
%       det_averaging: array specifying which detectors to average, dimension (2, number of averaged detectors)
%       average_time_segments: either an array specifying start and end timepoints of each chunk of timecourse to average
%           in which case dimension would be (number of segments, 2)
%           or 0 flag, under which case averaging time course parameters will be taken into account
%       decimate_factor: average time course parameter, factor to decimate by
%       moving_mean_window_length: average time course parameter, moving mean window length
%       avg_span: average time course parameter, number of data points to fit
%   dcsdatastruct:
%       counts: array with photon counts, dimension (timepoint, number of detectors)
%       g2: array with g2 values, dimension (timepoint, number of detectors, tau)
%       tau: array with tau values, dimension (1, tau)
%       t: array with timepoints corresponding to timepoints in g2 and counts, dimension (1, number of timepoints)
%
% output:
%   g2_data: array with autocorrelations averaged with time and detectors
%       dimension (ntau,ndetectors,ntimepoints)
%   intensities: array with photon counts
%       dimension (1, ntimepoints)
%   time_arr: array with timepoints
%       dimension (1, ntimepoints)
%   tau: array with tau values
%       dimension (1, ntau)
%
% author: Melissa Wu, <mwu22@mgh.harvard.edu>
% this function is part of the mcgeometry toolbox,
%(https://github.com/wumelissa/mc_geometry)
%% setting variables and detector averaging

corrset=dcsdatastruct.g2;
raw_intensities=dcsdatastruct.counts;
temp_time_arr=dcsdatastruct.t;
tau=dcsdatastruct.tau;

% average the detectors
if isempty(dcs_file.det_averaging)
    corrset_detavgd=corrset;
else
    corrset_detavgd=nan(size(corrset,1),size(dcs_file.det_averaging,1),size(corrset,3));
    for ch_idx=1:length(dcs_file.det_averaging)
        corrset_detavgd(:,ch_idx,:)=mean(corrset(:,dcs_file.det_averaging{ch_idx},:),2);
        corrset_detavgd_intensities(:,ch_idx)=mean(raw_intensities(:,dcs_file.det_averaging{ch_idx}),2);
    end
end

%% time segment averaging

% if user wants to average time segments (e.g. of baseline/intervention)
if isfield(dcs_file,'average_time_segments') && length(dcs_file.average_time_segments)>1
    
    % converting seconds to indices
    index_segment_array=nan(size(dcs_file.average_time_segments));
    
    for row=1:size(index_segment_array,1)
        for column=1:size(index_segment_array,2)
            [~,index_segment_array(row,column)]=min(abs(temp_time_arr-dcs_file.average_time_segments(row,column)));
        end
    end
    
    % averaging timepoints together
    for segment=1:size(dcs_file.average_time_segments,1)
        data_to_fit(segment,:,:)=mean(corrset_detavgd(index_segment_array(segment,1):index_segment_array(segment,2),:,:),1);
        time_arr(segment)=mean(temp_time_arr(index_segment_array(segment,1):index_segment_array(segment,2)));
        intensities(segment,:)=mean(corrset_detavgd_intensities(index_segment_array(segment,1):index_segment_array(segment,2),:),1);
    end
    
elseif isfield(dcs_file,'average_time_segments') && dcs_file.average_time_segments==0
    
    % creating final data array to fit - downsamples or averages the data as needed
    
    % decimating data
    for ch=1:size(corrset_detavgd,2)
        temp_intensity_arr=squeeze(raw_intensities(:,ch));
        decimate_intensity_arr(:,ch)=decimate(temp_intensity_arr,dcs_file.decimate_factor);
        for tau_point=1:size(corrset_detavgd,3)
            temp_arr=squeeze(corrset_detavgd(:,ch,tau_point));
            decimate_temp_arr=decimate(temp_arr,dcs_file.decimate_factor);
            decimate_corrset(:,ch,tau_point)=decimate_temp_arr;
        end
    end
    decimate_times=decimate(temp_time_arr,dcs_file.decimate_factor);
    
    % making a moving mean
    movmean_data=movmean(decimate_corrset,dcs_file.moving_mean_window_length,1);
    movmean_intensities=movmean(decimate_intensity_arr,dcs_file.moving_mean_window_length,1);
    movmean_times=decimate_times;
    
    % averaging timepoints together
    avg_data_to_fit=dcs_file.avg_span;
    last_time_seg=floor(size(movmean_data,1)/avg_data_to_fit)*avg_data_to_fit;
    time_arr(1)=0;
    idx=1;
    for I=1:avg_data_to_fit:last_time_seg
        data_to_fit(idx,:,:)=mean(movmean_data(I:(I+avg_data_to_fit-1),:,:),1);
        time_arr(idx)=mean(movmean_times(I:(I+avg_data_to_fit-1)));
        intensities(idx,:)=mean(movmean_intensities(I:(I+avg_data_to_fit-1),:),1);
        idx=idx+1;
    end
else
    error('Please enter a valid input option for dcs_file.average_time_segments\n');
end

%%

g2_data=permute(data_to_fit,[3 2 1]);

difft=diff(time_arr);
subtext=['Min sample duration: ' num2str(min(difft),'%.2f') ' s, max sample duration: ' num2str(max(difft),'%.2f') ' s'];

print_box('DCS DATA PRE-PROCESSED',subtext,100);
