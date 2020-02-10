function [BFi_arr,beta_arr]=analytical_fit_dcs(g2_data,tau,fit_options)

% analytical_fit_dcs returns fitted BFi and beta from an array of g2 (autocorrelation) curves, taus, and input parameter set
% 
% [BFi_arr,beta_arr]=analytical_fit_dcs(g2_data,tau,fit_options)
%
% input:
%   g2_data: array of autocorrelation values for each tau, detector, and timepoint, dimension (ntau,ndetect,ntimepoints)
%   tau: array of tau values for the autocorrelation, dimension (1,ntau)
%   fit_options: structure with fields:
%       x0: array with values [beta_initial_guess,bfi_initial_guess*1e9], dimension (1,2)
%       rhos_arr: array with values [distance1 distance2 ...] in mm, dimension (1,number of detectors)
%       lb: array with lower bounds for beta and BFi*1e9, dimension (1,2)
%       ub: array with upper bounds for beta and BFi*1e9, dimension (1,2)
%       ft: index of starting tau value of g2 curve to fit, dimension (1,1)
%       lt: index of starting tau value of g2 curve to fit, dimension (1,1)
%       lsq_options: output of optimoptions, specifying optimization options for lsqcurvefit fn in MATLAB
%       debug_plot: 0 or 1 flag for watching the fit plotted on a figure, dimension (1,1)
%
% output:
%   BFi_arr: array of BFi values for each detector channel and timepoint, dimension (ntimepoints,ndetect)
%   beta_arr: array of beta values for each detector channel and timepoint,
%   dimension (ntimepoints,ndetect)

% author: Melissa Wu, <mwu22@mgh.harvard.edu>
% this function is part of the mcgeometry toolbox,
%(https://github.com/wumelissa/mc_geometry)
%
%% setting variables

x0=fit_options.x0; 
rhos=fit_options.rhos_arr;
lb=fit_options.lb;
ub=fit_options.ub;

ft=fit_options.ft;
lt=fit_options.lt;
options=fit_options.lsq_options;

%% fitting for BFi 

text1='FITTING BFI ANALYTICALLY';
divider_width=100;
side1=ceil((divider_width-length(text1))/2);
side2=floor((divider_width-length(text1))/2);
fprintf([repmat('.',1,side1) text1 repmat('.',1,side2) '\n']);

reverseStr='';
idx=1;
for chidx=1:size(g2_data,2)
    for tpt_idx=1:size(g2_data,3)
        fit_options.rho=rhos(chidx);
        curr_g2=squeeze(g2_data(:,chidx,tpt_idx));
        test_x = lsqcurvefit(@(x,taus)semi_infinite_g2(x,tau(ft:lt),fit_options),x0,tau(ft:lt),curr_g2(ft:lt),lb,ub,options);
        if fit_options.debug_plot,
            computed_g2=semi_infinite_g2(test_x,tau(ft:lt),fit_options);
            figure(100); 
            semilogx(tau,curr_g2,tau(ft:lt),computed_g2);
            title([chidx tpt_idx]); 
            pause(.2);
        end
        beta_arr(tpt_idx,chidx)=test_x(1);
        BFi_arr(tpt_idx,chidx)=test_x(2)/1e9; % returns BFi in mm^2/s
        
        percentDone=100*idx/(size(g2_data,2)*size(g2_data,3));
        msg=sprintf('Percent done: %3.1f\n',percentDone);
        fprintf([reverseStr,msg])
        reverseStr=repmat(sprintf('\b'),1,length(msg));
        idx=idx+1;
    end
end