function plot_mc_fitting_result_with_tcd(BFi_arr,beta_arr,time_arr,analytical_fit_options,tcd_struct,save_plot_fullname);
%
% plot_mc_fitting_result_with_tcd(BFi_arr,beta_arr,time_arr,analytical_fit_options,tcd_struct,save_plot_fullname);
%
% plot Monte Carlo fitting results
% 
% input:
%   BFi_arr: array of fitted BFi values
%       either dimension (nchannels, ntimepoints)
%       or dimension (n superficial thicknesses, n middle thicknesses, nchannels, ntimepoints)
%   beta_arr: array of fitted beta values
%       either dimension (nchannels, ntimepoints)
%       or dimension (n superficial thicknesses, n middle thicknesses, nchannels, ntimepoints)
%   rhos_arr: array with source-detector distances in mm, dimension (1, nchannels)
%   region_splits: cell array denoting all combinations of tissue thicknesses to concatenate - see output of function
%       get_region_splits.m
%   save_plot_fullname: full filename to save plot
%   tcd_struct: structure with TCD data
%       t_left: left TCD timepoints
%       t_right: right TCD timepoints
%       left_min: left TCD lower envelope
%       left_max: left TCD upper envelope
%       right_min: right TCD lower envelope
%       right_max: right TCD upper envelope

% author: Melissa Wu, <mwu22@mgh.harvard.edu>
% this function is part of the mcgeometry toolbox,
%(https://github.com/wumelissa/mc_geometry)
%%

baseline_period=analytical_fit_options.baseline_period;
rhos_arr=analytical_fit_options.rhos_arr;

for det=1:length(rhos_arr)
    leg_arr{det}=[num2str(rhos_arr(det)) ' mm'];
end

%%

relative_tcd_time=tcd_struct.t_left-tcd_struct.t_left(1);

for idx=1:length(baseline_period)
    [~,baseline_period_indices_optical(idx)]=min(abs(time_arr-baseline_period(idx)));
    [~,baseline_period_tcd(idx)]=min(abs(relative_tcd_time-baseline_period(idx)));
end

normalized_left_min=tcd_struct.left_min/mean(tcd_struct.left_min(baseline_period_tcd(1):baseline_period_tcd(2)));
normalized_right_min=tcd_struct.right_min/mean(tcd_struct.right_min(baseline_period_tcd(1):baseline_period_tcd(2)));
normalized_left_max=tcd_struct.left_max/mean(tcd_struct.left_max(baseline_period_tcd(1):baseline_period_tcd(2)));
normalized_right_max=tcd_struct.right_max/mean(tcd_struct.right_max(baseline_period_tcd(1):baseline_period_tcd(2)));

%%

if length(size(BFi_arr))==2
    normalized_optical=BFi_arr./repmat(mean(BFi_arr(:,baseline_period_indices_optical(1):baseline_period_indices_optical(2)),2),[1 size(BFi_arr,2)]);
    
    fig1=figure;
    subplot(211)
    plot(time_arr,BFi_arr)
    grid on
    xlabel('time (seconds)'); ylabel('BF_i mm^2/s')
    legend('scalp BF_i','cerebral BF_i');
    title('BFi')
    
    subplot(212)
    plot(time_arr,beta_arr)
    grid on
    xlabel('time (seconds)'); ylabel('beta')
    legend(leg_arr);
    title('beta')
    drawnow
    
    fig2=figure;
    plot(time_arr,movmean(normalized_optical',3,1),'LineWidth',2);
    hold on
    plot(tcd_struct.t_left,movmean(normalized_left_min,3,1))
    plot(tcd_struct.t_right,movmean(normalized_right_min,3,1))
    plot(tcd_struct.t_left,movmean(normalized_left_max,3,1))
    plot(tcd_struct.t_right,movmean(normalized_right_max,3,1))
    grid on
    xlabel('seconds'); ylabel('rBFi')
    title('relative BFi')
    legend({'scalp BFi','cerebral BFi','normalized left TCD lower envelope','normalized left TCD upper envelope',...
    'normalized right TCD lower envelope','normalized right TCD upper envelope'})
    drawnow
    
    if ~isempty(save_plot_fullname)
        savefig(fig1,[save_plot_fullname '_mc_BFi_timecourse.fig']);
        savefig(fig2,[save_plot_fullname '_mc_beta_timecourse.fig']);
    end
else
    all_lengths=cellfun(@length,region_splits);
    sup_thickness_arr=unique(all_lengths(:,1));
    mid_thickness_arr=unique(all_lengths(:,2));
    
    [x,y]=find_subplot_dims(length(mid_thickness_arr));
    
    normalized_optical=BFi_arr./repmat(mean(BFi_arr(:,:,:,baseline_period_indices_optical(1):baseline_period_indices_optical(2)),4),[1 1 1 size(BFi_arr,4)]);
    
    for sup_thickness=sup_thickness_arr'
        h=figure;
        g=figure;
        
        idx=1;
        for mid_thickness=mid_thickness_arr'
            
            figure(h);
            subplot(x,y,idx)
            plot(time_arr,squeeze(BFi_arr(sup_thickness,mid_thickness,:,:))');
            grid on
            xlabel('seconds'); ylabel('BF_i mm^2/s');
            title(['mid thickness ' num2str(mid_thickness) ' mm']);
            drawnow
            
            single_optical_array=squeeze(normalized_optical(sup_thickness,mid_thickness,:,:));
            figure(g);
            subplot(x,y,idx)
            plot(time_arr,movmean(single_optical_array,3,2))
            hold on
            plot(tcd_struct.t_left,movmean(normalized_left_min,3,1))
            plot(tcd_struct.t_right,movmean(normalized_right_min,3,1))
            plot(tcd_struct.t_left,movmean(normalized_left_max,3,1))
            plot(tcd_struct.t_right,movmean(normalized_right_max,3,1))
            grid on
            xlabel('seconds'); ylabel('rBFi')
            title('relative BFi')
            drawnow
            
            idx=idx+1;
        end
        
        figure(h);
        subplot(x,y,1)
        legend('scalp BF_i','cerebral BF_i');
        suptitle(['BFi: superficial thickness ' num2str(sup_thickness) ' mm']);
        
        figure(g);
        subplot(x,y,1)
        legend({'scalp BFi','cerebral BFi', 'normalized left TCD lower envelope','normalized left TCD upper envelope',...
            'normalized right TCD lower envelope','normalized right TCD upper envelope'})
        suptitle(['rBFi with moving mean: superficial thickness ' num2str(sup_thickness) ' mm']);
    end
end