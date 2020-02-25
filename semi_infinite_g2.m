function y = semi_infinite_g2(x,taus,fit_options)
% returns autocorrelation using semi-infinite diffusion correlation
% equation
%
% input:
%   x: array with [beta, Db] value
%   taus: array of taus
%   fit_options: structure with fields
%       mu_a: absorption, mm-1
%       mu_s: scattering, mm-1
%       alpha: probability of scattering from a moving scatterer
%       lambda_dcs: wavelength, mm-1
%       rho: source-detector distance, mm

% author: Melissa Wu, <mwu22@mgh.harvard.edu>
% this function is part of the mcgeometry toolbox,
%(https://github.com/wumelissa/mc_geometry)

%% constant values, variable assignment

mu_a = fit_options.mu_a;
mu_s = fit_options.mu_s;
alpha = fit_options.alpha;
lambda_dcs = fit_options.lambda_dcs;
rho = fit_options.rho;
beta = x(1);
D = x(2)/1e9;

%%
%equations

n0=fit_options.n;
R = -1.440./n0^2+0.710/n0+0.668+0.0636*n0; % Effective reflection coefficient
ze = 2/3*(1+R)/(1-R); 
z_b = ze/mu_s;
z_0 = (mu_a + mu_s)^-1;
r_1 =(rho^2 + z_0^2)^0.5;
r_2 = (rho^2 + (z_0 + 2*z_b)^2)^0.5;
k_0 = 2*pi*fit_options.n/lambda_dcs;
delta_r2 = 6*D*taus;
K = (3*mu_a*mu_s + alpha*mu_s^2*k_0^2*delta_r2).^(0.5);
K_tau0 = (3*mu_a*mu_s)^0.5;
G1 = (3*mu_s/(4*pi))*(exp(-K*r_1)/r_1 - exp(-K*r_2)/r_2);
G1_tau0=(3*mu_s/(4*pi))*(exp(-K_tau0*r_1)/r_1 - exp(-K_tau0*r_2)/r_2);
small_g1 = abs(G1./G1_tau0);
y = (1+beta*(small_g1.^2))'; %final g2 equation

end

