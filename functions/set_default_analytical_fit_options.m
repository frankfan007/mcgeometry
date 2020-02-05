function analytical_fit_options=set_default_analytical_fit_options(baseline_period,save_plot);

% -------------------------------------------------------------------------
% EDIT: settings for BFi and beta fitting
% -------------------------------------------------------------------------

analytical_fit_options.baseline_period=baseline_period; % in seconds
analytical_fit_options.save_plot=save_plot;

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