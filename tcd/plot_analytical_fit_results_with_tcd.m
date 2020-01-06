function plot_analytical_fit_results_with_tcd(analytical_BFi,analytical_beta,intensities,time_arr,analytical_fit_options,tcd_struct)
%
% plot_analytical_fit_results_with_tcd(analytical_BFi,analytical_beta,intensities,time_arr,analytical_fit_options,tcd_struct)
%
% plot analytical fitting results

% input:
%   analytical_BFi: array with analytical BFi results, dimension (ntimepoints,nchannels)
%   analytical_beta: array with analytical beta results, dimension (ntimepoints,nchannels)
%   intensities: array with photon counts per second, dimension (ntimepoints,nchannels)
%   time_arr: array with times, dimension (1, ntimepoints)
%   analytical_fit_options: structure with fields
%       rhos_arr: array with source-detector distances in mm, dimension (1, nchannels)
%       save_plot: 0 or 1 flag to save plot
%       save_filename: full filename to save figure
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

normalized_optical=analytical_BFi./repmat(mean(analytical_BFi(baseline_period_indices_optical(1):baseline_period_indices_optical(2),:),1),[size(analytical_BFi,1) 1]);
normalized_left_min=tcd_struct.left_min/mean(tcd_struct.left_min(baseline_period_tcd(1):baseline_period_tcd(2)));
normalized_right_min=tcd_struct.right_min/mean(tcd_struct.right_min(baseline_period_tcd(1):baseline_period_tcd(2)));
normalized_left_max=tcd_struct.left_max/mean(tcd_struct.left_max(baseline_period_tcd(1):baseline_period_tcd(2)));
normalized_right_max=tcd_struct.right_max/mean(tcd_struct.right_max(baseline_period_tcd(1):baseline_period_tcd(2)));

%%

% plot the results
h=figure;
subplot(211)
plot(time_arr,analytical_BFi)
legend(leg_arr)
grid on
xlabel('seconds'); ylabel('BFi mm^2/s')
title('analytical fit BFi')

subplot(212)
plot(time_arr,movmean(normalized_optical,3,1),'LineWidth',2);
hold on
plot(tcd_struct.t_left,movmean(normalized_left_min,3,1))
plot(tcd_struct.t_right,movmean(normalized_right_min,3,1))
plot(tcd_struct.t_left,movmean(normalized_left_max,3,1))
plot(tcd_struct.t_right,movmean(normalized_right_max,3,1))
grid on
xlabel('seconds'); ylabel('rBFi')
title('relative BFi')
legend([leg_arr 'normalized left TCD lower envelope' 'normalized left TCD upper envelope'...
    'normalized right TCD lower envelope' 'normalized right TCD upper envelope'])
drawnow

g=figure;
subplot(211)
plot(time_arr,analytical_beta)
legend(leg_arr)
grid on
xlabel('seconds'); ylabel('beta')
title('analytical fit beta')

subplot(212)
plot(time_arr,intensities)
legend(leg_arr)
grid on
xlabel('seconds'); ylabel('photon counts/second')
title('intensities')
drawnow

%%

if analytical_fit_options.save_plot
    savefig(h,[analytical_fit_options.save_filename '_analytical_BFi_timecourse.fig']);
    savefig(g,[analytical_fit_options.save_filename '_analytical_beta_intensity.fig']);
end