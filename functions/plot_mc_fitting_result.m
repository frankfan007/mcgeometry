function plot_mc_fitting_result(BFi_arr,beta_arr,time_arr,rhos_arr,region_splits,save_plot_fullname);

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

%%

for det=1:length(rhos_arr)
    leg_arr{det}=[num2str(rhos_arr(det)) ' mm'];
end

if length(size(BFi_arr))==2
    fig1=figure;
    plot(time_arr,BFi_arr)
    grid on
    xlabel('time (seconds)'); ylabel('BF_i mm^2/s')
    legend('scalp BF_i','cerebral BF_i');
    title('BFi')
    drawnow
    
    fig2=figure;
    plot(time_arr,beta_arr)
    grid on
    xlabel('time (seconds)'); ylabel('beta')
    legend(leg_arr);
    title('beta')
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
    
    for sup_thickness=sup_thickness_arr'
        h=figure;
        idx=1;
        for mid_thickness=mid_thickness_arr'
            
            figure(h);
            subplot(x,y,idx)
            plot(time_arr,squeeze(BFi_arr(sup_thickness,mid_thickness,:,:))');
            grid on
            xlabel('time (seconds)'); ylabel('BF_i mm^2/s');
            title(['mid thickness ' num2str(mid_thickness) ' mm']);
            
            idx=idx+1;
        end
        
        figure(h);
        subplot(x,y,1)
        legend('scalp BF_i','cerebral BF_i');
        suptitle(['BFi: superficial thickness ' num2str(sup_thickness) ' mm']);
        drawnow
    end
end