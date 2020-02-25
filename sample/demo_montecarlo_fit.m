% -------------------------------------------------------------------------
% fit DCS data with Monte Carlo forward simulation
% -------------------------------------------------------------------------

% volume related

% IF ONLY FITTING A SINGLE COMBINATION OF LAYER CONCATENATIONS
% comment this out and use sup_thickness_arr/mid_thickness_arr variables if fitting multiple
% concatenate_tissue_layers_array={1,2,3,4:5}; % EDIT
% concatenate_tissue_layers_array={1:6,7:14,15:22};

% IF FITTING MULTIPLE LAYER CONCATENATIONS
% comment this out and use concatenate_tissue_layers_arr variable if fitting single one
sup_thickness_arr=1:7;
mid_thickness_arr=4:13;

% fitting related
fit_options.hold_superficial=1; % flag to set superficial BFi to short separation analytical BFi
fit_options.tau_range=4:110;
fit_options.mu_a=[0.0075 0.0248 0.0287];
fit_options.indices_of_mc_dets_to_use=[1 2];
fit_options.indices_of_dcs_dets_to_use=[1 2];
fit_options.beta_initial_guess=[0.5 0.5];
fit_options.bfi_initial_guess=[8e-6]; 
fit_options.hold_layer_indices=[1 2]; % MUST BE CONSISTENT WITH HOLD_SUPERFICIAL FLAG
fit_options.hold_layer_values=[nan 2e-8]; % if hold_superficial flag is 1, put a nan in the first index
fit_options.lsq_options=optimoptions('lsqcurvefit','TolFun',1e-5,'TolX',1e-3,'MaxFunEvals',1e5,'MaxIter',1e5,'Display','off');

% upper and lower bounds for fitting
fit_options.beta_lb=[0 0];
fit_options.bfi_lb=[0];
fit_options.beta_ub=[0.6 0.6];
fit_options.bfi_ub=[inf];

fit_options.tpts_to_fit=[]; % set to empty if you want to fit all timepoints in g2 array; otherwise, specify timepoints to fit

% measurement related
fit_options.n=1.37;
fit_options.wave=850e-6;
fit_options.k0=2*pi*fit_options.n/fit_options.wave;

% processing related
fit_options.do_plot=0;
fit_options.use_gpu=1; % number of gpu to use
fit_options.save_plot=1;

%%

% read MC photon history file

for det=1:length(analytical_fit_options.rhos_arr);detector_names{det}=[num2str(analytical_fit_options.rhos_arr(det)) ' mm'];end
[his_data,photon_indices,photon_fractions_retained,num_tissue_layers]=load_history_file(dir_struct.input_filename,mc_param.max_saved_photons,detector_names);
mc_his.his_array=his_data;

if exist('concatenate_tissue_layers_array','var')
    region_splits=concatenate_tissue_layers_array;
elseif exist('sup_thickness_arr','var') && exist('mid_thickness_arr','var')
    region_splits=get_region_splits(sup_thickness_arr,mid_thickness_arr,num_tissue_layers);
end

if isempty(fit_options.tpts_to_fit)
    fit_options.tpts_to_fit=1:size(g2_data,3);
end

save_plot_fullname=[dir_struct.mc_save_dir filesep save_filename];

% fit
[BFi_arr,beta_arr,rmse_arr,output_stats_arr]=loop_and_fit_dcs_parallel(g2_data,fit_options,mc_his,region_splits,analytical_BFi);
plot_mc_fitting_result(BFi_arr,beta_arr,time_arr,rhos_arr,region_splits,analytical_BFi,analytical_beta,baseline_period,fit_options.save_plot,save_plot_fullname);

mc_his.photon_indices=photon_indices;
mc_his.num_tissue_layers=num_tissue_layers;