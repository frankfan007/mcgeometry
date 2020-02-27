% -------------------------------------------------------------------------
% fit DCS data with Monte Carlo forward simulation
% -------------------------------------------------------------------------

cd ..
addpath(genpath(pwd))

mcx_flag=1; % EDIT: 1 if MCX was used, 0 if tMCimg was used
mcx_path='C:\Users\wu_me\OneDrive\Documents\mcx\'; % path where MCX is stored, if needed

addpath(genpath(mcx_path))

load sample_data_analytical_fit.mat
%%

% reading history file
input_filename='LargeSlab_MultiLyr1mm_085_mus.inp';
rhos_arr=[5 30];
max_saved_photons=1e5;

% only if using tMCimg - number of tissue layers used in the input volume
num_tissue_layers=21;

%%

% tissue layer indices to concatenate
sup_thickness_arr=1:7;
mid_thickness_arr=4:13;

% fitting related
fit_options.tau=tau';
fit_options.hold_superficial=1; % flag to set superficial BFi to short separation analytical BFi
fit_options.tau_range=10:116;
fit_options.mu_a=[0.0075 0.0248 0.0287];
fit_options.indices_of_mc_dets_to_use=[1 2]; % indices of Monte Carlo sim detectors to use
fit_options.indices_of_dcs_dets_to_use=[1 2]; % indices of analytical detectors to use
fit_options.beta_initial_guess=[0.5 0.5]; % initial guess for beta
fit_options.bfi_initial_guess=[8e-6]; % initial guess for BFi
fit_options.hold_layer_indices=[1 2]; % which tissue layers to hold at certain BFi values; must be consistent with hold_superficial flag
fit_options.hold_layer_values=[nan 2e-8]; % if hold_superficial flag is 1, put a nan in the first index
fit_options.lsq_options=optimoptions('lsqcurvefit','TolFun',1e-5,'TolX',1e-3,'MaxFunEvals',1e5,'MaxIter',1e5,'Display','off');

% upper and lower bounds for fitting
fit_options.beta_lb=[0 0]; % one for each detector
fit_options.bfi_lb=[0]; % one for each tissue layer
fit_options.beta_ub=[0.6 0.6];
fit_options.bfi_ub=[inf];

fit_options.tpts_to_fit=[]; % set to empty if you want to fit all timepoints in g2 array; otherwise, specify timepoints to fit

% measurement related
fit_options.n=1.37; % index of refraction
fit_options.wave=850e-6; % wavelength
fit_options.k0=2*pi*fit_options.n/fit_options.wave;

% processing related
fit_options.do_plot=0; % plot the fit
fit_options.use_gpu=0; % number of gpu to use
fit_options.save_plot=1;

save_plot_fullname='sample_data_slab_mc_fit';
%%

% read MC photon history file

for det=1:length(rhos_arr)
    detector_names{det}=[num2str(rhos_arr(det)) ' mm'];
end

if mcx_flag
    [his_data,photon_indices,photon_fractions_retained,num_tissue_layers]=load_history_file(input_filename,max_saved_photons,detector_names);
else
    history_filename=[input_filename(1:end-3) filesep '.his'];
    [his_data,photon_indices,photon_fractions_retained]=read_single_history_file(history_filename,num_tissue_layers,max_photons);
end

mc_his.his_array=his_data;
mc_his.photon_indices=photon_indices;
mc_his.num_tissue_layers=num_tissue_layers;

if exist('concatenate_tissue_layers_array','var')
    region_splits=concatenate_tissue_layers_array;
elseif exist('sup_thickness_arr','var') && exist('mid_thickness_arr','var')
    region_splits=get_region_splits(sup_thickness_arr,mid_thickness_arr,num_tissue_layers);
end

if isempty(fit_options.tpts_to_fit)
    fit_options.tpts_to_fit=1:size(g2_data,3);
end

% fit with Monte Carlo

[BFi_arr,beta_arr,rmse_arr,output_stats_arr]=loop_and_fit_dcs_parallel(g2_data,fit_options,mc_his,region_splits,analytical_BFi);

save(['.' filesep 'results' filesep 'sample_data_slab_multi_layer_mc_fit.mat'],'BFi_arr','beta_arr',...
    'rmse_arr','output_stats_arr','mc_his','fit_options')
%%

baseline_period=1:30;

plot_mc_fitting_result(BFi_arr,beta_arr,time_arr,rhos_arr,region_splits,analytical_BFi,analytical_beta,baseline_period,fit_options.save_plot,save_plot_fullname);
