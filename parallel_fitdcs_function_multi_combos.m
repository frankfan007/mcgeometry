function [BFi_arr,beta_arr,rmse_arr,output_stats_arr]=parallel_fitdcs_function_multi_combos(layer_combo,region_splits,input_his_array,mc_his,fit_options,g2_data,analytical_BFi)
% quick function to parallelize fitting for multiple superficial/deep layer
% combos
%
% needs fine-tuning

%%

for layer=1:size(region_splits,2)
    concatenate_tissue_layers_array{layer}=region_splits{layer_combo,layer};
end

his_temp=concatenate_layers(input_his_array,concatenate_tissue_layers_array);
mc_his.his_array=his_temp;

if fit_options.hold_superficial
    BFi_arr=nan(length(fit_options.bfi_initial_guess)+1,length(fit_options.tpts_to_fit));
else
    BFi_arr=nan(length(fit_options.bfi_initial_guess),length(fit_options.tpts_to_fit));
end

loop_tpts=tic;
for tpt_idx=1:length(fit_options.tpts_to_fit)
    tpt=fit_options.tpts_to_fit(tpt_idx);
    g2=squeeze(g2_data(:,:,tpt));
    if fit_options.hold_superficial
        fit_options.hold_layer_values(1)=analytical_BFi(tpt,1);
        [db_fit,beta_fit,rmse,output]=fit_g2_MC(g2,mc_his,fit_options);
        BFi_arr(1,tpt_idx)=analytical_BFi(tpt,1);
        BFi_arr(2:end,tpt_idx)=db_fit;
    else
        [db_fit,beta_fit,rmse,output]=fit_g2_MC(g2,mc_his,fit_options);
        BFi_arr(:,tpt_idx)=db_fit;
    end
    beta_arr(:,tpt_idx)=beta_fit;
    rmse_arr(tpt_idx)=rmse;
    output_stats_arr(tpt_idx)=output;
end
fprintf('Layer combination %d out of %d: finished in %.2f sec\n',layer_combo,size(region_splits,1),toc(loop_tpts))