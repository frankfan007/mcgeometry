% -------------------------------------------------------------------------
% fitting DCS data using semi-infinite analytical model
% -------------------------------------------------------------------------

% DCS file information
dcs_file.dir=['..' filesep 'raw'];
dcs_file.filename='sample_data.mat'; % EDIT

% averaging g2, detectors, and time segments
dcs_file.g2freq=10;
dcs_file.det_averaging={1:1,2:4}; 
dcs_file.average_time_segments=0; % either set zero or give span of seconds you want to average

% average data information
% the following are parameters for averaging the timecourse, and will be processed in the following order:
% decimate first, then moving mean, then averaging
dcs_file.decimate_factor=1; 
dcs_file.moving_mean_window_length=1; 
dcs_file.avg_span=1; 

%% analytical fit options

% more likely to be changed
analytical_fit_options.mu_a = 0.01; % mm-1
analytical_fit_options.mu_s = 0.85; % mm-1
analytical_fit_options.lambda_dcs = 850*1e-6; % mm-1
analytical_fit_options.rhos_arr=[5 30]; % mm
analytical_fit_options.ft=10; % first tau
analytical_fit_options.lt=116; % last tau
analytical_fit_options.debug_plot=0; 

% less likely to be changed
analytical_fit_options.alpha = 1;
analytical_fit_options.n=1.37; % refractive index
analytical_fit_options.beta_initial_guess=0.5;
analytical_fit_options.bfi_initial_guess=2e-6; % mm^2/s
analytical_fit_options.x0=[analytical_fit_options.beta_initial_guess analytical_fit_options.bfi_initial_guess*1e9]; % beta, then Db
analytical_fit_options.lb = zeros(size(analytical_fit_options.x0)); % lower bound for fitting
analytical_fit_options.ub=[]; % upper bound for fitting
analytical_fit_options.lsq_options = optimoptions('lsqcurvefit','TolFun',1e-8,'Display','off'); % options for lsqcurvefit

%% preparing DCS data and fitting

% measurement filename full path
dcs_file.measurement_file=[dcs_file.dir filesep dcs_file.filename];

% prepare DCS data
dcsdatastruct=prepare_sampledcs_data(dcs_file.measurement_file);

% pre-processing DCS data
[g2_data,time_arr,intensities,tau]=pre_process_dcs_data(dcs_file,dcsdatastruct);

% fit for BFi and beta and plot
[analytical_BFi,analytical_beta]=analytical_fit_dcs(g2_data,tau,analytical_fit_options);

%% plotting parameters

rhos_arr=analytical_fit_options.rhos_arr;

for det=1:length(rhos_arr)
    leg_arr{det}=[num2str(rhos_arr(det)) ' mm'];
end

%% plot the results

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
