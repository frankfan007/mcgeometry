function plot_mc_fitting_result_with_asl(BFi_arr,beta_arr,time_arr,rhos_arr,region_splits,asl_file,baseline_period,save_plot_fullname);

% plot_mc_fitting_result_with_asl(BFi_arr,beta_arr,time_arr,rhos_arr,region_splits,asl_file,baseline_period,save_plot_fullname)
%
% plot Monte Carlo fitting results
% 
% author: Melissa Wu, <mwu22@mgh.harvard.edu>
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
%   asl_file: name of asl file to load
%
% this function is part of the mcgeometry toolbox,
% (https://github.com/wumelissa/mc_geometry)
%%

asl_structure=load(asl_file);

for det=1:length(rhos_arr)
    leg_arr{det}=[num2str(rhos_arr(det)) ' mm'];
end

%%

for idx=1:length(baseline_period)
    [~,baseline_period_indices_optical(idx)]=min(abs(time_arr-baseline_period(idx)));
    [~,baseline_period_indices_asl(idx)]=min(abs(asl_structure.time_array-baseline_period(idx)));
end

normalized_asl=asl_structure.all_avgd_asl/mean(asl_structure.all_avgd_asl(baseline_period_indices_asl(1):baseline_period_indices_asl(2)));

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
    plot(time_arr,movmean(normalized_optical,3,2))
    hold on
    plot(asl_structure.time_array,movmean(normalized_asl,3,1),'LineWidth',2)
    grid on
    xlabel('seconds'); ylabel('rBFi')
    title('relative BFi, moving mean taken')
    legend('superficial BFi','deep BFi','ASL')
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
            plot(asl_structure.time_array,movmean(normalized_asl,3,1),'LineWidth',2)
            grid on
            xlabel('seconds'); ylabel('rBFi')
            title(['mid thickness ' num2str(mid_thickness) ' mm'])
            drawnow
            
            idx=idx+1;
        end
        
        figure(h);
        subplot(x,y,1)
        legend('scalp BF_i','cerebral BF_i');
        suptitle(['BFi: superficial thickness ' num2str(sup_thickness) ' mm']);
        
        figure(g);
        subplot(x,y,1)
        legend('scalp BF_i','cerebral BF_i','ASL');
        suptitle(['rBFi with moving mean: superficial thickness ' num2str(sup_thickness) ' mm']);
    end
end