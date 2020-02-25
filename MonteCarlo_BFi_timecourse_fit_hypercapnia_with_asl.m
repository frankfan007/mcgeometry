%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % ====================================== MONTE CARLO BASED BFI FITTING SCRIPT ====================================== %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CURRENTLY EDITING: NO LONGER IN USE. SEE SAMPLES FOLDER FOR DEMOS OF
% PARTICULAR FUNCTIONS

%%

% This script fits diffuse correlation spectroscopy data against a Monte Carlo based photon migration forward simulation 
% The volume used in the Monte Carlo simulation can be chosen from either a multi-layer slab, a multi-layer sample head, 
% or a subject-specific MRI scan provided by the user. The multi-layer sample head was derived from a FreeSurfer sample 
% T1 structural scan
% and subsequently post-processed using iterative image erosion into a 22-layer head volume
% For the slab and sample head, each layer is 1 mm thick - adjustable parameters (among others) for the MC fwd simulation include
% optical properties for each tissue layer and source-detector locations on head surface

% Please run this section by section to edit parameters as needed

% author: Melissa Wu, <mwu22@mgh.harvard.edu>
% this function is part of the mcgeometry toolbox,
%(https://github.com/wumelissa/mc_geometry)

%%

close all hidden
clearvars

%% ================================================= PATH SETTINGS ================================================== %%

% -------------------------------------------------------------------------
% EDIT: path settings
% -------------------------------------------------------------------------

subject_id='013';
session='2';

% ================ freesurfer directory ================ %

dir_struct.freesurfer_dir='/autofs/cluster/freesurfer/centos6_x86_64/stable6_0_0/'; % this is only usable in Linux

% ============ raw data and working folders ============ %

% set basepaths
if isunix
    raw_data_serverpath='';
    wd_serverpath=''; % working directory server path
    iso2mesh_serverpath='';
    mcx_serverpath='';
    dcs_mc_toolbox_serverpath='';
elseif ispc
    raw_data_serverpath='\\blinky.nmr.mgh.harvard.edu';
    wd_serverpath='\\blinky.nmr.mgh.harvard.edu';
    iso2mesh_serverpath='\\blinky.nmr.mgh.harvard.edu';
    mcx_serverpath='\\thm.nmr.mgh.harvard.edu';
    dcs_mc_toolbox_serverpath='\\thm.nmr.mgh.harvard.edu';
else
    error('Operating system not supported, please use a Linux or Windows system\n');
end

% working folder where processed data and volume will be saved
working_folder_unix=['/space/blinky/4/users/hypercapnia/Subject' subject_id '/Session' session '/'];
working_folder=[wd_serverpath working_folder_unix]; % hack for Windows
if ~exist(working_folder,'dir'), mkdir(working_folder);end

% raw data folder where DCS and MRI data is stored
raw_data_basepath=[raw_data_serverpath '/space/blinky/4/users/hypercapnia/'];
dir_struct.dcs_file_dir=[raw_data_basepath 'Subject' subject_id '/Session' session '/CW-DCS/']; % directory where DCS data is stored
dir_struct.mr_dir=[raw_data_basepath 'Subject' subject_id '_recon/mri/']; % directory where T1 is stored
dir_struct.asl_dir=[raw_data_basepath 'Subject' subject_id filesep 'perfusion' filesep];

% subject-based MRI volume directory, created in working folder
dir_struct.volume_dir=[working_folder filesep 'volume'];
dir_struct.volume_dir_unix=[working_folder_unix filesep 'volume'];
if ~exist(dir_struct.volume_dir), mkdir(dir_struct.volume_dir);end

% folder to save subject-specific Monte Carlo forward simulation 
dir_struct.mc_save_dir=[working_folder filesep 'mc'];
if ~exist(dir_struct.mc_save_dir,'dir'), mkdir(dir_struct.mc_save_dir);end

% ============= software and toolbox paths ============= %

% path where MCX software is stored
mcx_path='/space/thm/4/users/dibbyan/mcx/';
dir_struct.mcx_bin=[mcx_path 'bin/mcx'];

% path where iso2mesh software is stored
iso2mesh_path=[iso2mesh_serverpath '/space/blinky/4/users/hypercapnia/iso2mesh-2018-allinone/iso2mesh/']; 

% path where Monte Carlo geometry toolbox is stored
dir_struct.dcs_mc_toolbox_unix='/space/thm/4/users/melissa/mcgeometry/';
dir_struct.dcs_mc_toolbox=[dcs_mc_toolbox_serverpath dir_struct.dcs_mc_toolbox_unix];

% ==================== adding paths ==================== %

addpath(genpath(iso2mesh_path));
addpath([mcx_serverpath mcx_path '/utils']);
addpath(genpath(dir_struct.dcs_mc_toolbox));

% ============ full timecourse or heat plot ============ %

% don't edit
gen_options.full_timecourse=1;
gen_options.heat_plot=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% =============================================== DCS DATA SETTINGS ================================================ %%

% -------------------------------------------------------------------------
% EDIT: settings for DCS file and preparing g2s
% -------------------------------------------------------------------------

% DCS file information
dcs_file.dcsraw=0;
dcs_file.fastdcs=1;
dcs_file.filename='press_mod_1_run003_g2.mat'; % EDIT
dcs_file.g2freq=1;
dcs_file.det_averaging={1:1,2:4}; 

% average data information
% the following are parameters for averaging the timecourse, and will be processed in the following order:
% decimate first, then moving mean, then averaging
dcs_file.decimate_factor=1; 
dcs_file.moving_mean_window_length=1; 
dcs_file.avg_span=10; 

% -------------------------------------------------------------------------
% EDIT: settings for BFi and beta fitting
% -------------------------------------------------------------------------

% for normalizing and comparing to ASL
analytical_fit_options.baseline_period=[1 60]; % in seconds
analytical_fit_options.save_plot=1;

% more likely to be changed
analytical_fit_options.mu_a = 0.01; % mm-1
analytical_fit_options.mu_s = 1; % mm-1
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ==================================== MONTE CARLO FORWARD SIMULATION SETTINGS ===================================== %%

% -------------------------------------------------------------------------
% EDIT: specify volume, input file, and MC forward simulation parameters
% -------------------------------------------------------------------------

% which volume to use
volume_cfg.multi_layer_slab=0;
volume_cfg.multi_layer_head=1;
volume_cfg.subj_specific_mri=0;

% input and volume file names
mc_param.inp_filename='stairhead_5_30_mm_085_rotate_0'; % REQUIRED; write without extension .inp ; will use input file if exists or create one under variable name
mc_param.volume_name_noext=''; % volume name without extension - leave EMPTY if not using subject-specific MRI volume

% for wrapping probe around volume - only needed if creating new input file
ref_param.has_fiducial=0; % for subject-specific MRI generation, leave 0 if using multi-layer head or slab
ref_param.use_default_fiducial=1; % flag: 1 to use default fiducial location and 0 to choose your own
ref_param.det_distances=[5 30]; % mm
ref_param.fiducial_pos=[]; % if not empty, will automatically use this value as fiducial position
ref_param.rotate_deg=0; % if not empty, will automatically use this value as rotational degree

% for input file generation
mc_param.det_radii=[0.7 1]; % in mm
mc_param.mus=0.85; % either one value for each tissue layer, or single value for all tissue layers
mc_param.mua=0.01; % single value for all tissue layers
mc_param.n=1.37; 
mc_param.g=0.01; 
mc_param.num_phot_launched=1e9; 

% for Monte Carlo forward simulation - please see MCX documentation
mc_param.gpu_number=1; 
mc_param.max_detected_photons=3e6; % max number of detected photons 

% maximum number of saved photons from the MC history file
mc_param.max_saved_photons=1e5; % number of detected photons to save

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% =============================================== ASL FILE SETTINGS ================================================ %%

asl_filename='005_SURROUND.mat';

full_asl_filename=[dir_struct.asl_dir filesep asl_filename];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ========================================== MONTE CARLO BFI FIT SETTINGS ========================================== %%

% -------------------------------------------------------------------------
% EDIT: settings for Monte Carlo fitting
% -------------------------------------------------------------------------

% file related
[~,dcs_filename,dcs_fileext]=fileparts(dcs_file.filename);
save_filename=[dcs_filename '_' mc_param.inp_filename];

% volume related

% IF ONLY FITTING A SINGLE COMBINATION OF LAYER CONCATENATIONS
% comment this out and use sup_thickness_arr/mid_thickness_arr variables if fitting multiple
% concatenate_tissue_layers_array={1,2,3,4:5}; % EDIT
% concatenate_tissue_layers_array={1:5,6:10,11:22};

% IF FITTING MULTIPLE LAYER CONCATENATIONS
% comment this out and use concatenate_tissue_layers_arr variable if fitting single one
sup_thickness_arr=1:7;
mid_thickness_arr=4:13;

% fitting related
fit_options.hold_superficial=1; % flag to set superficial BFi to short separation analytical BFi
fit_options.tau_range=10:115;
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
fit_options.use_gpu=1;
fit_options.save_plot=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ============================================ PROCESSING - DO NOT EDIT ============================================ %%

% -------------------------------------------------------------------------
% preparing DCS data, analytical fitting, and plotting final BFi, beta, and intensities
% -------------------------------------------------------------------------

% measurement filename full path
dcs_file.measurement_file=[dir_struct.dcs_file_dir filesep dcs_file.filename];

% extra settings for DCS data
if gen_options.heat_plot + gen_options.full_timecourse ~=1
    error('Please specify either heat plot or full timecourse\n');
end
if gen_options.full_timecourse
    dcs_file.average_time_segments=0;
end

% EDIT: PROJECT SPECIFIC DCS DATA PREPARATION FOR DCS PRE-PROCESSING
if dcs_file.fastdcs
    dcsdatastruct=prepare_fastdcs_data(dcs_file);
elseif dcs_file.dcsraw
    dcsdatastruct=prepare_dcsraw_data(dcs_file);
end

% pre-processing DCS data
[g2_data,time_arr,intensities,tau]=pre_process_dcs_data(dcs_file,dcsdatastruct);

dcs_file.g2_data=g2_data;
dcs_file.time_arr=time_arr;
dcs_file.intensities=intensities;
dcs_file.tau=tau';
fit_options.tau=tau';

% fit for BFi and beta and plot
[analytical_BFi,analytical_beta]=analytical_fit_dcs(g2_data,tau,analytical_fit_options);

analytical_fit_options.save_filename=[dir_struct.dcs_file_dir filesep dcs_filename];
plot_analytical_fit_results_with_asl(analytical_BFi,analytical_beta,intensities,time_arr,analytical_fit_options,full_asl_filename)

% -------------------------------------------------------------------------
% read or generate input file, wraps probe, runs, and reads MC sim
% -------------------------------------------------------------------------

% setting some volume parameters
[dir_struct,ref_param,mc_param]=set_parameters_for_given_volume(dir_struct,ref_param,mc_param,volume_cfg);

% find input file status
input_file_status=return_input_file_status(dir_struct.input_filename);

% generates input file and runs MC simulation as needed
if input_file_status==-1
    if volume_cfg.subj_specific_mri
        generate_subject_seg(dir_struct,ref_param,mc_param)
    end
    volume_direc=dir_struct.volume_dir;
    toolbox_direc=[dir_struct.dcs_mc_toolbox filesep 'volumes'];
    if ref_param.has_fiducial
        volume_structure=generate_mri_mesh(volume_direc,...
            toolbox_direc,mc_param.volume_name,mc_param.meshvol_name); % load volume and generate mesh
    else
        volume_structure=generate_mri_mesh(volume_direc,...
            toolbox_direc,mc_param.volume_name,mc_param.meshvol_name); % load volume and generate mesh
    end
    ref_param=merge_structures(ref_param,volume_structure);
    ref_param=wrap_probe(ref_param); % wrap probe 
    generate_input_file(dir_struct,ref_param,mc_param); % generate input file
    run_mc_simulation(dir_struct,mc_param); % run MC forward simulation
elseif input_file_status==0
    run_mc_simulation(dir_struct,mc_param); % run MC forward simulation
elseif input_file_status==1
    subtext='reading history file';
    print_box('Input and history file exist',subtext,100);
end

% read MC photon history file
for det=1:length(analytical_fit_options.rhos_arr);detector_names{det}=[num2str(analytical_fit_options.rhos_arr(det)) ' mm'];end
[his_data,photon_indices,photon_fractions_retained,num_tissue_layers]=load_history_file(dir_struct.input_filename,mc_param.max_saved_photons,detector_names);
mc_his.his_array=his_data;
mc_his.photon_indices=photon_indices;
mc_his.num_tissue_layers=num_tissue_layers;

% -------------------------------------------------------------------------
% fit data against MC forward simulation
% -------------------------------------------------------------------------

% obtaining all combinations of superficial and middle thicknesses
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
plot_mc_fitting_result_with_asl(BFi_arr,beta_arr,time_arr,analytical_fit_options.rhos_arr,region_splits,full_asl_filename,analytical_BFi,analytical_beta,analytical_fit_options.baseline_period,fit_options.save_plot,save_plot_fullname);

save([save_plot_fullname '_full_timecourse.mat'],...
    'dir_struct','gen_options','full_asl_filename',...
    'analytical_fit_options','analytical_BFi','analytical_beta',...
    'BFi_arr','beta_arr','rmse_arr','time_arr',...
    'dcs_file','dcsdatastruct',...
    'volume_cfg',...
    'ref_param',...
    'mc_param','mc_his',...
    'fit_options','region_splits','BFi_arr','beta_arr','rmse_arr','output_stats_arr')
    