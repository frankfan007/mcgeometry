function plot_analytical_fit_results(analytical_BFi,analytical_beta,intensities,time_arr,analytical_fit_options)
% plot analytical fitting results
%
% input:
%   analytical_BFi: array with analytical BFi results, dimension (ntimepoints,nchannels)
%   analytical_beta: array with analytical beta results, dimension (ntimepoints,nchannels)
%   intensities: array with photon counts per second, dimension (ntimepoints,nchannels)
%   time_arr: array with times, dimension (1, ntimepoints)
%   analytical_fit_options: structure with fields
%       rhos_arr: array with source-detector distances in mm, dimension (1, nchannels)
%       save_plot: 0 or 1 flag to save plot
%       save_filename: full filename to save figure

% author: Melissa Wu, <mwu22@mgh.harvard.edu>
% this function is part of the mcgeometry toolbox,
%(https://github.com/wumelissa/mc_geometry)
%%

rhos_arr=analytical_fit_options.rhos_arr;

for det=1:length(rhos_arr)
    leg_arr{det}=[num2str(rhos_arr(det)) ' mm'];
end

% plot the results
h=figure;
plot(time_arr,analytical_BFi)
legend(leg_arr)
grid on
xlabel('seconds'); ylabel('BFi mm^2/s')
title('analytical fit BFi')
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