function [Dbfit,betafit,rmse,output_stats]=fit_g2_MC(g2_data,mc_his,fit_options)

%% fits for BFi and beta against data from a Monte Carlo forward simulation

% input:
%   g2_data: array with autocorrelations, dimension (ntau,number of detectors)
%   mc_his: structure with fields:
%       his_array: photon history array outputted from Monte Carlo simulation
%           array with path length and momentum transfer information stored for each detected photon
%           dimension (number of detected photons, 1 + 2*(number of tissue layers))
%           column 1 stores the index of the detector that the photon hit
%           columns 2:(1+number of tissue layers) stores the path lengths of the photons for each tissue layer
%           columns (2+number of tissue layers):end stores the momentum transfers of the photons for each tissue layer
%       photon_indices: cell with reference indices for his_array for photons detected by each detector
%           dimension (1, number of detectors)
%           ex: {1:1700,1701:2500,2501:3000}  denotes that columns (photons) 1:1700 in his_array are from detector 1,
%           columns (photons) 1701:2500 in his_array are from detector 2,
%           and columns (photons) 2501:3000 in his_array are from detector 3
%   fit_options: structure with fields:
%       tau_range: range (indices) of tau values for each curve to fit, dimension (1,length of range)
%       tau: array of tau values, dimension (1, ntau)
%       do_plot: 0 or 1 flag to plot fitted g2 curve on top of experimental data
%       hold_layer_indices: tissue layer indices of tissues that will have a constant BFi value (won't fit for BFi for those tissue layers)
%           dimension (1, number of tissues for which BFi is held)
%       hold_layer_values: constant BFi values for the tissue layers referenced in hold_layer_indices
%           dimension (1, number of tissues for which BFi is held)
%       mc_dets: detector indices of the detectors from the Monte Carlo forward simulation that will be used, dimension (1,number of detectors to be used)
%       dcs_dets: detector indices of the detectors from the DCS experimental data that will be used, dimension (1,number of detectors to be used)
%       beta_initial_guess: initial guess for beta, dimension (1, number of detectors)
%       bfi_initial_guess: initial guess for BFi, dimension (1, number of tissue layers that will be used)       

% output:
%   Dbfit: BFi for each tissue layer, dimension (1, number of tissue layers)
%   betafit: beta for each detector, dimension (1, number of detectors)
%   rmse: root mean square error (squared norm of the residual, where residual = fun(x,xdata) - ydata), dimension (1,1)
%   output_stats: output information about optimization process

%% defining and pre-allocating variables

his_array=mc_his.his_array;
photon_indices=mc_his.photon_indices;

tau_range=fit_options.tau_range;
tau=fit_options.tau;
do_plot=fit_options.do_plot;
hold_layer_indices=fit_options.hold_layer_indices;
hold_layer_values=fit_options.hold_layer_values*1e6;

num_layers=(size(his_array,2)-1)/2;
mc_dets=fit_options.indices_of_mc_dets_to_use;
dcs_dets=fit_options.indices_of_dcs_dets_to_use;
num_dets=length(mc_dets);

options=fit_options.lsq_options;

beta_x0=fit_options.beta_initial_guess;
db_x0=fit_options.bfi_initial_guess;
x0=[beta_x0 db_x0*1e6];

lb=[fit_options.beta_lb fit_options.bfi_lb];
ub=[fit_options.beta_ub fit_options.bfi_ub];

mts=cell(num_dets,1);
exps=cell(num_dets,1);
g1_norms=cell(num_dets,1);

%% setting GPU array (if applicable)

if fit_options.use_gpu~=0
    gpuDevice(fit_options.use_gpu);
    tau=transpose(gpuArray(tau(tau_range)));
    mu_a=gpuArray(fit_options.mu_a);
    k0=gpuArray(fit_options.k0);
elseif fit_options.use_gpu==0
    tau=transpose(tau(tau_range));
    mu_a=fit_options.mu_a;
    k0=fit_options.k0;
else
    error('fit_options.use_gpu value should be either 1 (use GPU) or 0 (use CPU)\n')
end

%% calculating momentum transfer and g1 for each detector

for detector=1:num_dets
    idx_start=photon_indices{mc_dets(detector)}(1);
    idx_end=photon_indices{mc_dets(detector)}(2);
    mtransfer=his_array(idx_start:idx_end,(num_layers+2):end);
    path_length=his_array(idx_start:idx_end,2:(num_layers+1));
    exp1=repmat(exp(-mu_a*path_length'),[length(tau) 1])';
    g1_norm=sum(exp(-mu_a*path_length'));
    mts{detector}=mtransfer;
    exps{detector}=exp1;
    g1_norms{detector}=g1_norm;
end

g2_fit=squeeze(g2_data(tau_range,dcs_dets));

%% fit

singlefit=tic;
[b_db_fit,rmse,~,~,output_stats]=lsqcurvefit(@(x,tau)mc_g2_Db_beta_hold_layer(x,tau,g2_fit,num_dets,mts,exps,g1_norms,num_layers,hold_layer_indices,hold_layer_values,k0,do_plot),x0,tau,g2_fit,lb,ub,options);
output_stats.time=toc(singlefit);

betafit=b_db_fit(1:num_dets);
Dbfit=b_db_fit((num_dets+1):end)/1e6;
