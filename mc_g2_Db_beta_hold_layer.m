function result=mc_g2_Db_beta_hold_layer(s,tau,expg2,num_dets,mts,exps,g1_norms,num_layers,hold_layer_indices,hold_layer_values,k0,varargin)
% returns g2 for each detector named in the chs variable
% must receive either single Db value, or one Db for each tissue type
% same for mu_a and mu_sp
% if more than one channel is supplied in expg1, they should be concatenated column
% vectors
%
% input:
%   s: array with 1:num_dets values as initial guesses for beta, and num_dets+1:end values as initial guesses or hold values for BFi
%       dimension (1, number of detectors + number of tissue layers)
%   tau: array with tau values
%   expg2: expected g2 to plot on top of result, dimension (1, ntau)
%   num_dets: number of detectors
%   mts: cell with array of momentum transfers for each photon for each tissue layer, dimension (1, number of detectors)
%   exps: cell array, exponential of absorption multiplied by path length, dimension (1, number of detectors)
%   g1_norms: cell array, sum of exponential of absorption multiplied by path length, dimension (1, number of detectors)
%   num_layers: number of tissue layers
%   hold_layer_indices: indices of tissue layers to hold BFi for, dimension (1, number of hold layers)
%   hold_layer_values: values of BFi to hold, dimension (1, number of hold layers)
%   k0: wavenumber at DCS wavelength
%   varargin: 
%       varargin{1}: flag to show fit figure        

% author: Melissa Wu, <mwu22@mgh.harvard.edu>
% this function is part of the mcgeometry toolbox,
%(https://github.com/wumelissa/mc_geometry)
%%

% parsing out the inputs and setting the Db "hold" values if there are any
beta=s(1:num_dets);
if isempty(hold_layer_values)
    Db=transpose(s((num_dets+1):end))/1e6;
else
    fit_Db=s(num_dets+1:end);
    fit_idx=1;
    hold_idx=1;
    for layer=1:num_layers
        if ismember(layer,hold_layer_indices)
            Db(layer)=hold_layer_values(hold_idx);
            hold_idx=hold_idx+1;
        else
            Db(layer)=fit_Db(fit_idx);
            fit_idx=fit_idx+1;
        end
    end
    Db=Db'/1e6;
end

%% adding up g2

if ~isempty(varargin), showfig=varargin{1}; else, showfig=0; end

if length(Db)==1, Db=ones(num_layers,1)*Db; end

for I=1:num_dets
    temp_big=exp(-(k0.^2.*2*mts{I}*(Db*tau))).*exps{I};
    g1(:,I)=sum(temp_big)'/g1_norms{I};
    g2(:,I)=1+beta(I)*(g1(:,I).^2);
end
  
g2=gather(g2);
result=g2;
 
 %% plotting
 
 if showfig
     figure(149);
     hold off;
     semilogx(repmat(tau,[1 num_dets]),expg2(:));
     hold on
     semilogx(repmat(tau,[1 num_dets]),g2(:),'r');
     title({sprintf('%0.2f %0.2f %0.2f\n',beta), sprintf('%1.2e %1.2e %1.2e %1.2e\n',Db([1 end]))})
     ylim([0.8 1.7]);
     drawnow
 end
 