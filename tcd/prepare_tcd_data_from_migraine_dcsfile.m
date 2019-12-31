function [tcd_struct,output_answer,varargout]=prepare_tcd_data_from_migraine_dcsfile(dcs_file,time_arr,varargin)

% specifically to process TCD data from Juliette and Phoebe's migraine data

% input
%   dcs_file: structure with fields
%       measurement_file: full filename of measurement file to load
%   time_arr: array with timepoints, dimension (1,ntimepoints)
%   varargin:
%       varargin{1}: array with flags to use left or right side TCD, dimension (1,2)

% output
%   tcd_struct: structure with TCD data
%       t_left: left TCD timepoints
%       t_right: right TCD timepoints
%       left_min: left TCD lower envelope
%       left_max: left TCD upper envelope
%       right_min: right TCD lower envelope
%       right_max: right TCD upper envelope
%   output_answer: array with flags to use left or right side TCD, dimension (1,2)
%   varargout:
%       varargout{1}: array denoting baseline and hypercapnia periods, dimension (nperiods,2)

%%

load(dcs_file.measurement_file)

%%

if exist('tcd','var') && ~isempty(tcd) && exist('t_nirs','var')
    if ~isempty(varargin)
        input_answer=varargin{1};
        [tcd_cell,output_answer]=get_tcd_envelope(t_nirs,tcd,input_answer);
    else
        [tcd_cell,output_answer]=get_tcd_envelope(t_nirs,tcd);
    end
    
    t_left_raw=tcd_cell{1};
    left_min_raw=tcd_cell{2};
    left_max_raw=tcd_cell{3};
    
    t_right_raw=tcd_cell{4};
    right_min_raw=tcd_cell{5};
    right_max_raw=tcd_cell{6};
    
    if isfield(dcs_file,'average_time_segments') && length(dcs_file.average_time_segments)>1
        
        left_index_segment_array=nan(size(dcs_file.average_time_segments));
        right_index_segment_array=nan(size(dcs_file.average_time_segments));
        
        if ~isempty(t_left_raw)
            for row=1:size(dcs_file.average_time_segments,1)
                for column=1:size(dcs_file.average_time_segments,2)
                    [~,left_index_segment_array(row,column)]=min(abs(t_left_raw-dcs_file.average_time_segments(row,column)));
                end
            end
            for segment=1:size(dcs_file.average_time_segments,1)
                left_min(segment)=median(left_min_raw(left_index_segment_array(segment,1):left_index_segment_array(segment,2)),'omitnan');
                left_max(segment)=median(left_max_raw(left_index_segment_array(segment,1):left_index_segment_array(segment,2)),'omitnan');
                t_left(segment)=median(t_left_raw(left_index_segment_array(segment,1):left_index_segment_array(segment,2)),'omitnan');
            end
        else
            left_min=[];
            left_max=[];
            t_left=[];
        end
        
        if ~isempty(t_right_raw)
            for row=1:size(dcs_file.average_time_segments,1)
                for column=1:size(dcs_file.average_time_segments,2)
                    [~,right_index_segment_array(row,column)]=min(abs(t_right_raw-dcs_file.average_time_segments(row,column)));
                end
            end
            for segment=1:size(dcs_file.average_time_segments,1)
                right_min(segment)=median(right_min_raw(right_index_segment_array(segment,1):right_index_segment_array(segment,2)),'omitnan');
                right_max(segment)=median(right_max_raw(right_index_segment_array(segment,1):right_index_segment_array(segment,2)),'omitnan');
                t_right(segment)=median(t_right_raw(right_index_segment_array(segment,1):right_index_segment_array(segment,2)),'omitnan');
            end
        else
            right_min=[];
            right_max=[];
            t_right=[];
        end
        varargout{1}=[];
        
    elseif isfield(dcs_file,'average_time_segments') && dcs_file.average_time_segments==1
        
        trig_indices=find(s_allHC>0);
        nb_trig=length(trig_indices);
        trig_seconds=round(t_nirs(trig_indices));
        
        if ~isempty(t_left_raw)
            segment_array=nan(2*length(trig_seconds),2);
            segment_array_in_seconds=nan(2*length(trig_seconds),2);
            for segment=1:2:2*length(trig_seconds)
                idx=floor(segment/2)+1;
                first_bound_val=trig_seconds(idx)-30;
                second_bound_val=trig_seconds(idx);
                third_bound_val=trig_seconds(idx)+30;
                
                [~,first_bound_index]=min(abs(t_left_raw-first_bound_val));
                [~,second_bound_index]=min(abs(t_left_raw-second_bound_val));
                [~,third_bound_index]=min(abs(t_left_raw-third_bound_val));
                
                segment_array(segment,1)=first_bound_index;
                segment_array(segment,2)=second_bound_index;
                segment_array(segment+1,1)=second_bound_index;
                segment_array(segment+1,2)=third_bound_index;
                
                segment_array_in_seconds(segment,1)=first_bound_val;
                segment_array_in_seconds(segment,2)=second_bound_val;
                segment_array_in_seconds(segment+1,1)=second_bound_val;
                segment_array_in_seconds(segment+1,2)=third_bound_val;
            end
            
            for segment=1:size(segment_array,1)
                left_min(segment)=median(left_min_raw(segment_array(segment,1):segment_array(segment,2)),'omitnan');
                left_max(segment)=median(left_max_raw(segment_array(segment,1):segment_array(segment,2)),'omitnan');
                t_left(segment)=median(t_left_raw(segment_array(segment,1):segment_array(segment,2)),'omitnan');
            end
        else
            left_min=[];
            left_max=[];
            t_left=[];
        end
        
        if ~isempty(t_right_raw)
            segment_array=nan(2*length(trig_seconds),2);
            for segment=1:2:2*length(trig_seconds)
                idx=floor(segment/2)+1;
                first_bound_val=trig_seconds(idx)-30;
                second_bound_val=trig_seconds(idx);
                third_bound_val=trig_seconds(idx)+30;
                
                [~,first_bound_index]=min(abs(t_right_raw-first_bound_val));
                [~,second_bound_index]=min(abs(t_right_raw-second_bound_val));
                [~,third_bound_index]=min(abs(t_right_raw-third_bound_val));
                
                segment_array(segment,1)=first_bound_index; %% THESE ARE CURRENTLY IN SECONDS, NOT THE ACTUAL INDICES!!!
                segment_array(segment,2)=second_bound_index;
                segment_array(segment+1,1)=second_bound_index;
                segment_array(segment+1,2)=third_bound_index;
                
                segment_array_in_seconds(segment,1)=first_bound_val;
                segment_array_in_seconds(segment,2)=second_bound_val;
                segment_array_in_seconds(segment+1,1)=second_bound_val;
                segment_array_in_seconds(segment+1,2)=third_bound_val;
            end
            
            for segment=1:size(segment_array,1)
                right_min(segment)=median(right_min_raw(segment_array(segment,1):segment_array(segment,2)),'omitnan');
                right_max(segment)=median(right_max_raw(segment_array(segment,1):segment_array(segment,2)),'omitnan');
                t_right(segment)=median(t_right_raw(segment_array(segment,1):segment_array(segment,2)),'omitnan');
            end
        else
            right_min=[];
            right_max=[];
            t_right=[];
        end
        varargout{1}=segment_array_in_seconds;
        
    elseif isfield(dcs_file,'average_time_segments') && dcs_file.average_time_segments==0

        avg_seconds_to_fit=median(diff(time_arr));
        
        if ~isempty(t_left_raw)
            decimate_envelope_data{1}=decimate(left_min_raw(~isnan(left_min_raw)),round(avg_seconds_to_fit/median(diff(t_left_raw))));
            decimate_envelope_data{2}=decimate(left_max_raw(~isnan(left_max_raw)),round(avg_seconds_to_fit/median(diff(t_left_raw))));
            decimate_envelope_data{3}=decimate(t_left_raw(~isnan(left_max_raw)),round(avg_seconds_to_fit/median(diff(t_left_raw))));
        else
            decimate_envelope_data{1}=[];
            decimate_envelope_data{2}=[];
            decimate_envelope_data{3}=[];
        end
        
        if ~isempty(t_right_raw)
            decimate_envelope_data{4}=decimate(right_min_raw(~isnan(right_min_raw)),round(avg_seconds_to_fit/median(diff(t_right_raw))));
            decimate_envelope_data{5}=decimate(right_max_raw(~isnan(right_max_raw)),round(avg_seconds_to_fit/median(diff(t_right_raw))));
            decimate_envelope_data{6}=decimate(t_right_raw(~isnan(right_min_raw)),round(avg_seconds_to_fit/median(diff(t_right_raw))));
        else
            decimate_envelope_data{4}=[];
            decimate_envelope_data{5}=[];
            decimate_envelope_data{6}=[];
        end
        
        t_left=decimate_envelope_data{3};
        t_right=decimate_envelope_data{6};
        left_min=decimate_envelope_data{1};
        left_max=decimate_envelope_data{2};
        right_min=decimate_envelope_data{4};
        right_max=decimate_envelope_data{5};
        
        varargout{1}=[];
    else
        error('Please enter a valid input option for dcs_file.average_time_segments');
    end
   
    tcd_struct.t_left=t_left;
    tcd_struct.t_right=t_right;
    tcd_struct.left_min=left_min;
    tcd_struct.left_max=left_max;
    tcd_struct.right_min=right_min;
    tcd_struct.right_max=right_max;
    
else
    error('No TCD data available\n');
end
