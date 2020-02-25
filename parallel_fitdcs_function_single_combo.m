function [BFi,beta_fit,rmse,output]=parallel_fitdcs_function_single_combo(mc_his,fit_options,g2_data,analytical_BFi,tpt)
% quick function to parallelize fitting for a single superficial/deep layer
% combos
%
% needs fine-tuning

%%

g2=squeeze(g2_data(:,:,tpt));

if fit_options.hold_superficial
    BFi=nan(length(fit_options.bfi_initial_guess)+1,1);
    fit_options.hold_layer_values(1)=analytical_BFi(tpt,1);
    [db_fit,beta_fit,rmse,output]=fit_g2_MC(g2,mc_his,fit_options);
    BFi(1)=analytical_BFi(tpt,1);
    BFi(2:end)=db_fit;
else
    BFi=nan(length(fit_options.bfi_initial_guess),1);
    [db_fit,beta_fit,rmse,output]=fit_g2_MC(g2,mc_his,fit_options);
    BFi=db_fit;
end