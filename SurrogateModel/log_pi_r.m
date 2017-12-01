%% Last update: 02/06/2017

% PROJECT NAME: Multi-variate UQ
% SUBPROJECT:   Multi-output emulator

% This function gives prior \pi_{\boldsymbol{R}}(\boldsymbol{r}) for the roughness parameters r_i
% It is a normalized scale of the product of i.i.d. vague (albeit proper) Log-Logistic prior type 
% LaTeX formula: $\pi_{\boldsymbol{R}}(\boldsymbol{r})=\prod_{i=1}^p (1+r_i^2)^{-1}$

function [ p ] = log_pi_r( r )

p = -sum(log(1+r.^2));

end


