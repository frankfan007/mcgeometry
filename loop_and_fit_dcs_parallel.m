function [BFi_arr,beta_arr,rmse_arr,output_stats_arr]=loop_and_fit_dcs_parallel(g2_data,fit_options,mc_his,region_splits,analytical_BFi)

% loops through timepoints of array with g2 data and fits each timepoint for BFi and beta
%
% input:
%   g2_data: array with autocorrelations, dimension (ntau,number of detectors,number of timepoints)
%   fit_options: substructure with fields
%       hold_superficial: 0 or 1 flag to hold superficial layer to short separation analytical fit BFi
%       bfi_initial_guess: initial guess for BFi, dimension (1, number of tissue layers)
%       tpts_to_fit: array of indices of timepoints to fit, dimension (1, timepoints to fit)
%       hold_layer_values: array of BFi values to hold tissue BFi at, if any are held
%           dimension (1, number of held layers)
%   mc_his: structure with fields
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
%   region_splits: cell array denoting all combinations of tissue thicknesses to concatenate - see output of function
%       get_region_splits.m
%       dimension (number of combinations of layers to split,1)
%   analytical_BFi: array with BFi calculated from semi-infinite analytical fit, dimension (number of timepoints,number of detectors)
%       
% output:
%   BFi_arr: array of fitted BFi values, dimensionality depends on if fitting for multiple combinations of superficial/middle
%       thicknesses, or if fitting for just one
%       if fitting for just one: dimension (number of detectors, number of timepoints)
%       if fitting for multiple: dimension (maximum superficial thickness, maximum middle thickness, number of tissue layers, number of timepoints)
%   beta_arr: array of fitted beta values, dimensionality depends on if fitting for multiple combinations of superficial/middle
%       thicknesses, or if fitting for just one
%       if fitting for just one: dimension (number of detectors, number of timepoints)
%       if fitting for multiple: dimension (maximum superficial thickness, maximum middle thickness, number of tissue layers, number of timepoints)
%   rmse_arr: array of root mean squared errors from the fit, dimensionality depends on if fitting for multiple combinations of superficial/middle
%       thicknesses, or if fitting for just one
%       if fitting for just one: dimension (1,number of timepoints)
%       if fitting for multiple: dimension (maximum superficial thickness, maximum middle thickness, number of timepoints)
%   output_stats_arr: structure array of output stats from the fit, dimensionality depends on if fitting for multiple combinations of superficial/middle
%       thicknesses, or if fitting for just one
%       if fitting for just one: dimension (1,number of timepoints)
%       if fitting for multiple: dimension (maximum superficial thickness, maximum middle thickness, number of timepoints)

% author: Melissa Wu, <mwu22@mgh.harvard.edu>
% this function is part of the mcgeometry toolbox,
%(https://github.com/wumelissa/mc_geometry)
%% setting variables

input_his_array=mc_his.his_array;

%% loop through timepoints and fit for BFi and beta

text1='FITTING BFI AGAINST MONTE CARLO FORWARD SIM';
divider_width=100;
side1=ceil((divider_width-length(text1))/2);
side2=floor((divider_width-length(text1))/2);
fprintf([repmat('.',1,side1) text1 repmat('.',1,side2) '\n']);

tpts_to_fit=fit_options.tpts_to_fit;

if size(region_splits,1)==1
    
    concatenate_tissue_layers_array=region_splits;
    his_temp=concatenate_layers(mc_his.his_array,concatenate_tissue_layers_array);
    mc_his.his_array=his_temp;
    if fit_options.hold_superficial
        BFi_arr=nan(length(fit_options.bfi_initial_guess)+1,length(fit_options.tpts_to_fit));
    else
        BFi_arr=nan(length(fit_options.bfi_initial_guess),length(fit_options.tpts_to_fit));
    end
    
    parfor idx=1:length(fit_options.tpts_to_fit)
        tpt=tpts_to_fit(idx);
        fprintf('Fitting file %d out of %d \n ',idx,length(fit_options.tpts_to_fit))
        tic
        [BFi_arr(:,idx),beta_arr(:,idx),rmse_arr(:,idx),output_stats_arr(:,idx)]=parallel_fitdcs_function_single_combo(mc_his,fit_options,g2_data,analytical_BFi,tpt);
        fprintf('Finished in %2.2f seconds', toc)
    end
    
elseif size(region_splits,1)>1
    
    all_thicknesses=cellfun('size',region_splits,2);
    
    if fit_options.hold_superficial
        BFi_arr=nan(max(all_thicknesses(:,1)),max(all_thicknesses(:,2)),length(fit_options.bfi_initial_guess)+1,length(fit_options.tpts_to_fit));
        temp_BFi_arr=nan(size(region_splits,1),length(fit_options.bfi_initial_guess)+1,length(fit_options.tpts_to_fit));
    else
        BFi_arr=nan(max(all_thicknesses(:,1)),max(all_thicknesses(:,2)),length(fit_options.bfi_initial_guess),length(fit_options.tpts_to_fit));
        temp_BFi_arr=nan(size(region_splits,1),length(fit_options.bfi_initial_guess),length(fit_options.tpts_to_fit));
    end
    
    regions=region_splits;
        
    parfor layer_combo=1:size(region_splits,1)
        
        [temp_BFi_arr(layer_combo,:,:),temp_beta_arr(layer_combo,:,:),temp_rmse_arr(layer_combo,:),temp_output_stats_arr(layer_combo,:)]=parallel_fitdcs_function_multi_combos(layer_combo,...
            region_splits,input_his_array,mc_his,fit_options,g2_data,analytical_BFi);
    end
    
    for layer_combo=1:size(region_splits,1)
        
        sup_layer_thickness=length(regions{layer_combo,1});
        mid_layer_thickness=length(regions{layer_combo,2});
        
        BFi_arr(sup_layer_thickness,mid_layer_thickness,:,:)=temp_BFi_arr(layer_combo,:,:);
        beta_arr(sup_layer_thickness,mid_layer_thickness,:,:)=temp_beta_arr(layer_combo,:,:);
        rmse_arr(sup_layer_thickness,mid_layer_thickness,:)=temp_rmse_arr(layer_combo,:);
        output_stats_arr(sup_layer_thickness,mid_layer_thickness,:)=temp_output_stats_arr(layer_combo,:);
        
    end
end

