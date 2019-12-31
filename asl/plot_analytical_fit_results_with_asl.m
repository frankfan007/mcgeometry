function plot_analytical_fit_results_with_asl(analytical_BFi,analytical_beta,intensities,time_arr,analytical_fit_options,asl_file)

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
%   asl_file: name of asl file to load

% this function is part of the mcgeometry toolbox

%%

baseline_period=analytical_fit_options.baseline_period;
rhos_arr=analytical_fit_options.rhos_arr;

asl_structure=load(asl_file);

for det=1:length(rhos_arr)
    leg_arr{det}=[num2str(rhos_arr(det)) ' mm'];
end

%%

for idx=1:length(baseline_period)
    [~,baseline_period_indices_optical(idx)]=min(abs(time_arr-baseline_period(idx)));
    [~,baseline_period_indices_asl(idx)]=min(abs(asl_structure.time_array-baseline_period(idx)));
end

normalized_optical=analytical_BFi./repmat(mean(analytical_BFi(baseline_period_indices_optical(1):baseline_period_indices_optical(2),:),1),[size(analytical_BFi,1) 1]);
normalized_asl=asl_structure.all_avgd_asl/mean(asl_structure.all_avgd_asl(baseline_period_indices_asl(1):baseline_period_indices_asl(2)));

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
plot(time_arr,movmean(normalized_optical,3,1));
hold on
plot(asl_structure.time_array,movmean(normalized_asl,3,1),'LineWidth',2);
grid on
xlabel('seconds'); ylabel('rBFi')
title('relative BFi, moving mean taken')
legend([leg_arr 'ASL'])
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
    savefig(h,[analytical_fit_options.save_filename '_analytical_BFi_with_asl.fig']);
    savefig(g,[analytical_fit_options.save_filename '_analytical_beta_intensity.fig']);
end